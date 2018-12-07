{-# LANGUAGE UnicodeSyntax, PatternGuards, RecordWildCards #-}
import qualified Network.Socket as Net
import qualified System.Environment
import qualified Request
import qualified RequestEval
import qualified Sys
import qualified Data.Map as Map
import qualified Network.BSD
import qualified Cxx.Show
import qualified IRC
import qualified Data.ByteString

import IRC (Command(..), Prefix(..))
import Control.Exception (bracketOnError)
import System.IO (hSetBinaryMode, hFlush, Handle, IOMode(..), stdout)
import Control.Monad (forever, when)
import Control.Arrow (first)
import Control.Monad.State (execStateT, lift, StateT, get)
import Control.Monad.Writer (execWriterT, tell)
import System.Console.GetOpt (OptDescr(..), ArgDescr(..), ArgOrder(..), getOpt, usageInfo)
import System.Locale.SetLocale (setLocale, Category(..))
import Text.Regex (Regex, subRegex, mkRegex, mkRegexWithOpts) -- Todo: Text.Regex truncates Char's >256. Get rid of it.
import Data.Char (isSpace, isPrint, isDigit)
import Data.List (isSuffixOf, isPrefixOf)
import Data.Map (Map)
import Data.SetOps
import Debug.Trace (trace)
import Util ((.), elemBy, caselessStringEq, maybeM, describe_new_output,
  orElse, full_evaluate, withResource, mapState',
  strip_utf8_bom, none, takeBack, replaceInfix, classify_diagnostic)
import Sys (rate_limiter)

import Prelude hiding ((.))
import Prelude.Unicode hiding ((∈))

data IrcBotConfig = IrcBotConfig
  { server :: Net.HostName, port :: Net.PortNumber
  , password :: Maybe String
  , max_response_length :: Int
  , chans :: [String], key_chans :: [(String, String)]
  , nick :: String, nick_pass :: Maybe String, alternate_nick :: String
  , also_respond_to :: [String]
  , allow_nickless_requests_in :: [String]
  , blacklist :: [String]
  , no_output_msg :: String
  , channel_response_prefix :: String
      -- A first occurrence of the string "nick" is replaced with the nick of the requester.
  , join_trigger :: Maybe IRC.Message
      -- Defaults to RPL_WELCOME. Can be set to NickServ/cloak confirmations and such.
  , censor :: [Regex]
  , rate_limit_messages, rate_limit_window :: Int
  , serve_private_requests :: Bool
  , clang_by_default :: Bool
  } deriving Read

instance Read Regex where
  readsPrec i s = first (\r → mkRegexWithOpts r True False) . readsPrec i s

data Opt = Help deriving Eq

optsDesc :: [OptDescr Opt]
optsDesc = [Option "h" ["help"] (NoArg Help) "Display this help and exit."]

help :: String
help = usageInfo "Usage: sudo geordi-irc [option]...\nOptions:" optsDesc ++ "\nSee README.xhtml for more information."

getArgs :: IO [Opt]
getArgs = do
  args ← System.Environment.getArgs
  case getOpt RequireOrder optsDesc args of
    (_, _, err:_) → fail $ init err
    (_, w:_, []) → fail $ "superfluous command line argument: " ++ w
    (opts, [], []) → return opts

do_censor :: IrcBotConfig → String → String
do_censor cfg s = foldr (\r t → subRegex r t "<censored>") s (censor cfg)

main :: IO ()
main = do
 setLocale LC_ALL (Just "")
 opts ← getArgs
 if Help ∈ opts then putStrLn help else do
  cfg@IrcBotConfig{..} ← read . getContents
  full_evaluate $ do_censor cfg "abc" -- So that any mkRegex failures occur before we start connecting.
  putStrLn $ "Connecting to " ++ server ++ ":" ++ show port
  withResource (connect server (fromIntegral port)) $ \h → do
   putStrLn "Connected"
   evalRequest ← RequestEval.evaluator
   limit_rate ← rate_limiter rate_limit_messages rate_limit_window
   let send m = limit_rate >> IRC.send h (IRC.Message Nothing m)
   maybeM password $ send . Pass
   send $ Nick nick
   send $ User nick 0 nick
   flip execStateT (∅) $ forever $ do
     l ← lift $ Data.ByteString.hGetLine h
     case IRC.decode l of
       Nothing → lift $ putStrLn "Malformed IRC message."
       Just m → do
         lift $ print m
         r ← on_msg evalRequest cfg (Data.ByteString.length l == 511) m
         lift $ mapM_ print r >> hFlush stdout >> mapM_ send r
   return ()

discarded_lines_description :: Int → String
discarded_lines_description s =
  " [+ " ++ show s ++ " discarded line" ++ (if s == 1 then "" else "s") ++ "]"

describe_lines :: [String] → String
describe_lines [] = ""
describe_lines (x:xs)
  | xs == [] || classify_diagnostic x == Just "error" = x
  | otherwise = x ++ discarded_lines_description (length xs)

data ChannelMemory = ChannelMemory
  { context :: Request.Context
  , last_outputs :: [String]
  , last_nonrequest :: String }
type ChannelMemoryMap = Map String ChannelMemory

emptyChannelMemory :: IrcBotConfig → ChannelMemory
emptyChannelMemory IrcBotConfig{..} = ChannelMemory
  { context = Request.Context Cxx.Show.noHighlighting clang_by_default []
  , last_outputs = []
  , last_nonrequest = "" }

is_request :: IrcBotConfig → Where → String → Maybe String
is_request IrcBotConfig{..} _ s
  | Just (n, r) ← Request.is_addressed_request s
  , any (caselessStringEq n) (nick : alternate_nick : also_respond_to)
    = Just r
is_request IrcBotConfig{..} (InChannel c) s
  | elemBy caselessStringEq c allow_nickless_requests_in
  , Just r ← Request.is_nickless_request s
    = Just r
is_request _ Private s = Just s
is_request _ _ _ = Nothing

type Reason = String
data Permission = Allow | Deny (Maybe Reason)
data Where = Private | InChannel String

request_allowed :: IrcBotConfig → String → Maybe IRC.UserName → Maybe IRC.ServerName → Where → Permission
request_allowed cfg _ _ _ Private | not (serve_private_requests cfg) =
  Deny $ Just "This bot does not serve private requests."
request_allowed cfg nickname _ _ _ | nickname ∈ blacklist cfg = Deny Nothing
request_allowed _ _ _ _ _ = Allow

type Eraser = String → Maybe String

digits :: Eraser
digits (x : y : s) | isDigit x, isDigit y = Just s
digits (x : s) | isDigit x = Just s
digits _ = Nothing

color_code :: Eraser
color_code ('\x3' : ',' : s) = digits s
color_code ('\x3' : s) = case digits s of
  Just (',' : s') → digits s'
  Just s' → Just s'
  Nothing → Just s
color_code _ = Nothing

apply_eraser :: Eraser → String → String
apply_eraser _ [] = []
apply_eraser p s@(h:t) = p s `orElse` (h : apply_eraser p t)

strip_color_codes :: String → String
strip_color_codes = apply_eraser color_code
  {- Todo: The above is *much* more naturally expressed as:
        subRegex r s "" where r = mkRegex "\x3(,[[:digit:]]{1,2}|[[:digit:]]{1,2}(,[[:digit:]]{1,2})?)?"
  Unfortunately, Text.Regex is broken: it truncates Char's, resulting in spurious matches. -}

version_response :: String
version_response = "Geordi C++ bot - http://www.eelis.net/geordi/"

strip_discord :: String → String
strip_discord s = do
    subRegex r s "" where r = mkRegex "^<[a-zA-Z0-9 ]+>"

on_msg :: (Functor m, Monad m) ⇒
  (String → Request.Context → [(String, String)] → m Request.Response) → IrcBotConfig → Bool → IRC.Message → StateT ChannelMemoryMap m [IRC.Command]
on_msg eval cfg@IrcBotConfig{..} full_size m@(IRC.Message prefix c) = execWriterT $ do
  when (join_trigger == Just m) join
  case c of
    Quit | Just (NickName n _ _) ← prefix, n == nick → send $ Nick nick
    PrivMsg _ "\1VERSION\1" | Just (NickName n _ _) ← prefix →
      send $ Notice n $ "\1VERSION " ++ version_response ++ "\1"
    NickNameInUse → send $ Nick alternate_nick
    Ping x → send $ Pong x
    PrivMsg _ ('\1':_) → return ()
    PrivMsg to txt' | Just (NickName who muser mserver) ← prefix → do
      let
        txt = filter isPrint $ strip_discord | trace("discord " ++ show txt) False = undefined $ strip_color_codes $ strip_utf8_bom txt'
        private = elemBy caselessStringEq to [nick, alternate_nick]
        w = if private then Private else InChannel to
        wher = if private then who else to
        reply s = send $ PrivMsg wher $ take max_response_length $
            (if private then id else (replaceInfix "nick" who channel_response_prefix ++)) $
            if null s then no_output_msg else do_censor cfg s
      mem@ChannelMemory{..} ← (`orElse` emptyChannelMemory cfg) . Map.lookup wher . lift get
      case (dropWhile isSpace . is_request cfg w txt) of
        Nothing → lift $ mapState' $ insert (wher, mem{last_nonrequest = txt'})
        Just r' → case request_allowed cfg who muser mserver w of
         Deny reason → maybeM reason reply
         Allow → do
          let r = if r' == "^" then last_nonrequest else r'
          if full_size ∧ none (`isSuffixOf` r) ["}", ";"] then reply $ "Request likely truncated after `" ++ takeBack 15 r ++ "`." else do
            -- The "}"/";" test gains a reduction in false positives at the cost of an increase in false negatives.
          let extra_env = [("GEORDI_REQUESTER", who), ("GEORDI_WHERE", wher)]
          Request.Response history_modification output ← lift $ lift $ eval r context extra_env
          let output' = describe_lines $ dropWhile null $ lines output
          let lo = take 50 last_outputs
          lift $ mapState' $ insert (wher, mem
            { context = maybe id Request.modify_history history_modification context
            , last_outputs = output' : lo })
          reply $ describe_new_output lo output'
    Welcome → do
      maybeM nick_pass $ send . PrivMsg "NickServ" . ("identify " ++)
      when (join_trigger == Nothing) join
    Invite _ _ → join
    _ → return ()
  where
    send = tell . (:[])
    join = send $ Join chans key_chans

connect :: Net.HostName → Net.PortNumber → IO Handle
  -- Mostly copied from Network.connectTo. We can't use that one because we want to set SO_KEEPALIVE (and related) options on the socket, which can't be done on a Handle.
connect host portn = do
 proto ← Network.BSD.getProtocolNumber "tcp"
 let hints = Net.defaultHints { Net.addrSocketType = Net.Stream, Net.addrProtocol = proto }
 target ← head . Net.getAddrInfo (Just hints) (Just host) (Just $ show portn)
 bracketOnError (Net.socket (Net.addrFamily target) Net.Stream proto) Net.sClose $ \sock → do
  Sys.setKeepAlive sock 30 10 5
  Net.connect sock (Net.addrAddress target)
  h ← Net.socketToHandle sock ReadWriteMode
  hSetBinaryMode h True
  return h

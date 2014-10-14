FROM ubuntu:14.04
RUN apt-get update && apt-get install -y libreadline-dev libboost-dev build-essential libgmp-dev pkg-config libseccomp-dev

ADD https://www.haskell.org/platform/download/2014.2.0.0/haskell-platform-2014.2.0.0-unknown-linux-x86_64.tar.gz /

RUN tar xvf haskell-platform-2014.2.0.0-unknown-linux-x86_64.tar.gz
RUN /usr/local/haskell/ghc-7.8.3-x86_64/bin/activate-hs
RUN rm /haskell-platform-2014.2.0.0-unknown-linux-x86_64.tar.gz
# 14.04's haskell-platform package is too old

RUN cabal update && cabal install --global filepath process deepseq mtl syb unix utf8-string network containers readline parsec Diff directory regex-compat base-unicode-symbols setops streams semigroups regex-posix template-haskell transformers pointed distributive comonad contravariant profunctors semigroupoids irc setlocale
# (We could rely on the list in geordi.cabal, but having these here shortens the development cycle when I'm testing changes in geordi.)

COPY etc /geordi/etc
COPY src /geordi/src

RUN cabal install --global /geordi/src --prefix=/usr

RUN /geordi/src/prep-image

CMD ["/usr/bin/geordi-local"]

WORKDIR /geordi/run
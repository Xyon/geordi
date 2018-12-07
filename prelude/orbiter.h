#include <limits>
#include "windef.h"
// ======================================================================
/// \defgroup constants Some useful general constants
// ======================================================================
//@{
const double PI    = 3.14159265358979323846;	///< pi
const double PI05  = 1.57079632679489661923;	///< pi/2
const double PI2   = 6.28318530717958647693;	///< pi*2
const double RAD   = PI/180.0;      ///< factor to map degrees to radians
const double DEG   = 180.0/PI;      ///< factor to map radians to degrees
const double C0    = 299792458.0;   ///< speed of light in vacuum [m/s]
const double TAUA  = 499.004783806; ///< light time for 1 AU [s]
const double AU    = C0*TAUA;       ///< astronomical unit (mean geocentric distance of the sun) [m]
const double GGRAV = 6.67259e-11;   ///< gravitational constant [m^3 kg^-1 s^-2]
const double G     = 9.81;          ///< gravitational acceleration [m/s^2] at Earth mean radius
const double ATMP  = 101.4e3;       ///< atmospheric pressure [Pa] at Earth sea level
const double ATMD  = 1.293;         ///< atmospheric density [kg/m^3] at Earth sea level
//@}

// ======================================================================
// API data types
// ======================================================================

class VESSEL;
class CELBODY;
class ExternMFD;
class Interpreter;

namespace oapi {
	class Module;
	class Sketchpad;
	class Font;
	class Pen;
	class Brush;

// ======================================================================
/// \defgroup defines Defines and Enumerations
/// \defgroup structures Structure definitions
// ======================================================================


// ======================================================================
/// \ingroup defines
/// \defgroup handle Handles
// ======================================================================
//@{
/// \brief Handle for objects (vessels, stations, planets)
typedef void *OBJHANDLE;

/// \brief Handle for visuals
typedef void *VISHANDLE;

/// \brief Handle for meshes
typedef void *MESHHANDLE;

/// \brief Handle for graphics-client-specific meshes
typedef int *DEVMESHHANDLE;
//struct DEVMESHHANDLE {
//	DEVMESHHANDLE() { hMesh = NULL; }
//	DEVMESHHANDLE(MESHHANDLE h) { hMesh = h; }
//	DWORD id;
//	MESHHANDLE hMesh;
//	operator int() { return (int)hMesh; }
//};

/// \brief Handle for bitmap surfaces and textures (panels and panel items)
typedef void *SURFHANDLE;

/// \brief Handle for 2D instrument panels
typedef void *PANELHANDLE;

/// \brief Handle for file streams
typedef void *FILEHANDLE;

/// \brief Handle for script interpreters
typedef void *INTERPRETERHANDLE;

/// \brief Handle for thrusters
typedef void *THRUSTER_HANDLE;

/// \brief Handle for logical thruster groups
typedef void *THGROUP_HANDLE;

/// \brief Propellant resource handle
typedef void *PROPELLANT_HANDLE;

/// \brief Handle for particle streams
typedef void *PSTREAM_HANDLE;

/// \brief Handle for vessel docking ports
typedef void *DOCKHANDLE;

/// \brief Handle vor vessel passive attachment points
typedef void *ATTACHMENTHANDLE;

/// \brief Handle for vessel airfoils
typedef void *AIRFOILHANDLE;

/// \brief Handle for vessel aerodynamic control surfaces
typedef void *CTRLSURFHANDLE;

/// \brief Handle for a navigation radio transmitter (VOR, ILS, IDS, XPDR)
typedef void *NAVHANDLE;

/// \brief Handle for animation components
typedef void *ANIMATIONCOMPONENT_HANDLE;

/// \brief Handle for custom items added to Launchpad "Extra" list
typedef void *LAUNCHPADITEM_HANDLE;

/// \brief Handle for onscreen annotation objects
typedef void *NOTEHANDLE;
//@}

typedef enum { FILE_IN, FILE_OUT, FILE_APP } FileAccessMode;
typedef enum { ROOT, CONFIG, SCENARIOS, TEXTURES, TEXTURES2, MESHES, MODULES } PathRoot;

// ===========================================================================
/**
 * \defgroup vec Vectors and matrices
 * Vectors and matrices are used to represent positions, velocities, translations,
 *   rotations, etc. in the 3-dimensional object space. Orbiter provides the
 *   %VECTOR3 and %MATRIX3 structures for 3-D vectors and matrices. A number
 *   of utility functions allow common operations such as matrix-vector
 *   products, dot and vector products, etc.
 */
// ===========================================================================
//@{
/**
 * \brief 3-element vector
 */
typedef union {
	double data[3];               ///< array data interface
	struct { double x, y, z; };   ///< named data interface
} VECTOR3;

/**
 * \brief 3x3-element matrix
 */
typedef union {
	double data[9];               ///< array data interface (row-sorted)
	struct { double m11, m12, m13, m21, m22, m23, m31, m32, m33; }; ///< named data interface
} MATRIX3;

typedef union {      // 4x4 matrix
	double data[16];
	struct { double m11, m12, m13, m14, m21, m22, m23, m24, m31, m32, m33, m34, m41, m42, m43, m44; };
} MATRIX4;
//@}

/**
 * \ingroup structures
 * \brief colour definition
 */
typedef struct {
	float r;    ///< read colour component [0..1]
	float g;    ///< green colour component [0..1]
	float b;    ///< blue colour component [0..1]
	float a;    ///< alpha (opacity) component (0..1)
} COLOUR4;

/**
 * \ingroup structures
 * \brief vertex definition including normals and texture coordinates
 */
typedef struct {
	float x;     ///< vertex x position
	float y;     ///< vertex y position
	float z;     ///< vertex z position
	float nx;    ///< vertex x normal
	float ny;    ///< vertex y normal
	float nz;    ///< vertex z normal
	float tu;    ///< vertex u texture coordinate
	float tv;    ///< vertex v texture coordinate
} NTVERTEX;

/**
 * \ingroup structures
 * \brief Defines a mesh group (subset of a mesh).
 *
 * A mesh group contains a vertex list, an index list,
 * a material and texture index, and a set of flags.
 */
typedef struct {
	NTVERTEX *Vtx;     ///< vertex list
	WORD *Idx;         ///< index list
	DWORD nVtx;        ///< vertex count
	DWORD nIdx;        ///< index count
	DWORD MtrlIdx;     ///< material index (>= 1, 0=none)
	DWORD TexIdx;      ///< texture index (>= 1, 0=none)
	DWORD UsrFlag;     ///< user-defined flag
	WORD zBias;        ///< z bias
	WORD Flags;        ///< internal flags
} MESHGROUP;

const DWORD MAXTEX = 1;  // max. extra textures per mesh group

/**
 * \ingroup structures
 * \brief extended mesh group definition
 */
typedef struct {
	NTVERTEX *Vtx;     ///< vertex list
	WORD *Idx;         ///< index list
	DWORD nVtx;        ///< vertex count
	DWORD nIdx;        ///< index count
	DWORD MtrlIdx;     ///< material index (>= 1, 0=none)
	DWORD TexIdx;      ///< texture index (>= 1, 0=none)
	DWORD UsrFlag;     ///< user-defined flag
	WORD zBias;        ///< z bias
	WORD Flags;        ///< internal flags
	DWORD TexIdxEx[MAXTEX]; ///< additional texture indices
	float TexMixEx[MAXTEX]; ///< texture mix values
} MESHGROUPEX;

/**
 * \ingroup defines
 * \defgroup grpedit Mesh group editing flags
 * These constants can be applied to the \e flags field of the
 *   \ref GROUPEDITSPEC structure to define which parts of a
 *   mesh group are to be modified.
 * \note The GRPEDIT_SETUSERFLAG, GRPEDIT_ADDUSERFLAG and
 *   GRPEDIT_DELUSERFLAG flags are mutually exclusive. Only one
 *   can be used at a time.
 * \sa GROUPEDITSPEC, oapiEditMeshGroup
 */
//@{
#define GRPEDIT_SETUSERFLAG 0x0001 ///< replace the group's UsrFlag entry with the value in the GROUPEDITSPEC structure.
#define GRPEDIT_ADDUSERFLAG 0x0002 ///< Add the UsrFlag value to the group's UsrFlag entry
#define GRPEDIT_DELUSERFLAG 0x0004 ///< Remove the UsrFlag value from the group's UsrFlag entry
#define GRPEDIT_VTXCRDX     0x0008 ///< Replace vertex x-coordinates
#define GRPEDIT_VTXCRDY     0x0010 ///< Replace vertex y-coordinates
#define GRPEDIT_VTXCRDZ     0x0020 ///< Replace vertex z-coordinates
#define GRPEDIT_VTXCRD      (GRPEDIT_VTXCRDX | GRPEDIT_VTXCRDY | GRPEDIT_VTXCRDZ) ///< Replace vertex coordinates
#define GRPEDIT_VTXNMLX     0x0040 ///< Replace vertex x-normals
#define GRPEDIT_VTXNMLY     0x0080 ///< Replace vertex y-normals
#define GRPEDIT_VTXNMLZ     0x0100 ///< Replace vertex z-normals
#define GRPEDIT_VTXNML      (GRPEDIT_VTXNMLX | GRPEDIT_VTXNMLY | GRPEDIT_VTXNMLZ) ///< Replace vertex normals
#define GRPEDIT_VTXTEXU     0x0200 ///< Replace vertex u-texture coordinates
#define GRPEDIT_VTXTEXV     0x0400 ///< Replace vertex v-texture coordinates
#define GRPEDIT_VTXTEX      (GRPEDIT_VTXTEXU | GRPEDIT_VTXTEXV) ///< Replace vertex texture coordinates
#define GRPEDIT_VTX         (GRPEDIT_VTXCRD | GRPEDIT_VTXNML | GRPEDIT_VTXTEX) ///< Replace vertices
//@}

/**
 * \ingroup structures
 * \brief Structure used by \ref oapiEditMeshGroup to define the
 *   group elements to be replaced.
 * \note Only the group elements specified in the \e flags entry will
 *   be replaced or modified. The elements that are to remain unchanged
 *   can be left undefined in the GROUPEDITSPEC structure. For example,
 *   if only GRPEDIT_VTXCRDX is specified, only the 'x' fields in the
 *   Vtx array need to be assigned.
 * \note to replace individual vertices in the group, the nVtx entry
 *   should contain the number of vertices to be replaced, the vIdx
 *   array should contain the indices (>= 0) of the vertices to be
 *   replaced, and Vtx should contain the new vertex values of those
 *   vertices. If vIdx==NULL, vertices are replaced in sequence from
 *   the beginning of the group's vertex list.
 *   nVtx must be less or equal the number of vertices in the group.
 * \sa oapiEditMeshGroup, grpedit
 */
typedef struct {
	DWORD flags;   ///< flags (see \ref grpedit)
	DWORD UsrFlag; ///< Replacement for group UsrFlag entry
	NTVERTEX *Vtx; ///< Replacement for group vertices
	DWORD nVtx;    ///< Number of vertices to be replaced
	WORD *vIdx;    ///< Index list for vertices to be replaced
} GROUPEDITSPEC;

/**
 * \ingroup structures
 * \brief material definition
 */
typedef struct {
	COLOUR4 diffuse;   ///< diffuse component
	COLOUR4 ambient;   ///< ambient component
	COLOUR4 specular;  ///< specular component
	COLOUR4 emissive;  ///< emissive component
	float power;       ///< specular power
} MATERIAL;

/**
 * \brief Kepler orbital elements
 *
 * A set of 6 scalar parameters defining the state of an object in a 2-body
 * (Keplerian) orbit. The orbital trajectory is a conic section, either
 * closed (circular, elliptic), or open (parabolic, hyperbolic).
 * \note semi-major axis a is positive for closed orbits, and negative for
 *   open orbits (in that case, a is referred to as real semi-axis).
 * \note eccentricity e:
 *   - circular orbit: e = 0
 *   - elliptic orbit: 0 < e < 1
 *   - parabolic orbit: e = 1
 *   - hyperbolic orbit: e > 1
 * \note The a and e parameters define the shape of the orbit, the i, theta
 *   and omegab parameters define the orientation of the orbital plane in
 *   space, and the L parameter defines the object position along the
 *   trajectory at a given time.
 * \note This is a generic data format. Additional data are required to
 *   fully define an object's state in space (position and velocity
 *   vectors). These include the position of the orbited body, the
 *   orientation of the reference coordinate system, and the date to which
 *   the mean longitude parameter refers.
 * \sa ORBITPARAM, \subpage orbit
 */
typedef struct {
	double a;          ///< semi-major axis [m]
	double e;          ///< eccentricity
	double i;          ///< inclination [rad]
	double theta;      ///< longitude of ascending node [rad]
	double omegab;     ///< longitude of periapsis [rad]
	double L;          ///< mean longitude at epoch
} ELEMENTS;

/**
 * \brief Secondary orbital parameters derived from the primary \ref ELEMENTS.
 *
 * This members of this structure provide additional parameters to the
 * primary elements of contained in the ELEMENTS structure.
 * \note SMi: for open orbits, this represents the imaginary semi-axis
 * \note PeD: distance to lowest point of the orbit from focal point
 * \note ApD: distance of highest point of the orbit from focal point. Only
 *   defined for closed orbits.
 * \note T: orbit period only defined for closed orbits.
 * \note PeT: For open orbits, this is negative after periapis passage
 * \note ApT: Only defined for closed orbits.
 * \sa ELEMENTS, \subpage orbit
 */
typedef struct {
	double SMi;        ///< semi-minor axis [m]
	double PeD;        ///< periapsis distance [m]
	double ApD;        ///< apoapsis distance [m]
	double MnA;        ///< mean anomaly [rad]
	double TrA;        ///< true anomaly [rad]
	double MnL;        ///< mean longitude [rad]
	double TrL;        ///< true longitude [rad]
	double EcA;        ///< eccentric anomaly [rad]
	double Lec;        ///< linear eccentricity [m]
	double T;          ///< orbit period [s]
	double PeT;        ///< time to next periapsis passage [s]
	double ApT;        ///< time to next apoapsis passage [s]
} ORBITPARAM;

/**
 * \ingroup structures
 * \brief Planetary atmospheric constants structure
 */
typedef struct {
	double p0;         ///<     pressure at mean radius ('sea level') [Pa]
	double rho0;       ///<     density at mean radius
	double R;          ///<     specific gas constant [J/(K kg)]
	double gamma;      ///<     ratio of specific heats, c_p/c_v
	double C;          ///<     exponent for pressure equation (temporary)
	double O2pp;       ///<     partial pressure of oxygen
	double altlimit;   ///<     atmosphere altitude limit [m]
	double radlimit;   ///<     radius limit (altlimit + mean radius)
	double horizonalt; ///<     horizon rendering altitude
	VECTOR3 color0;    ///<     sky colour at sea level during daytime
} ATMCONST;

/** \brief Atmospheric parameters structure */
typedef struct {
	double T;          ///<     temperature [K]
	double p;          ///<     pressure [Pa]
	double rho;        ///<     density [kg/m^3]
} ATMPARAM;

/** \brief Engine status */
typedef struct {
	double main;       ///<     -1 (full retro) .. +1 (full main)
	double hover;      ///<     0 .. +1 (full hover)
	int attmode;       ///<     0=rotation, 1=translation
} ENGINESTATUS;

/**
 * \ingroup defines
 * \defgroup exhaustflag Bitflags for EXHAUSTSPEC flags field.
 * \sa EXHAUSTSPEC
 */
//@{
#define EXHAUST_CONSTANTLEVEL 0x0001 ///< exhaust level is constant
#define EXHAUST_CONSTANTPOS   0x0002 ///< exhaust position is constant
#define EXHAUST_CONSTANTDIR   0x0004 ///< exhaust direction is constant
//@}

/**
 * \brief Engine exhaust render parameters
 * \sa VESSEL::AddExhaust(EXHAUSTSPEC*)
 */
typedef struct {
	THRUSTER_HANDLE th;///<  handle of associated thruster (or NULL if none)
	double *level;     ///<  pointer to variable containing exhaust level (0..1)
	VECTOR3 *lpos;     ///<  pointer to exhaust position vector [m]
	VECTOR3 *ldir;     ///<  pointer to engine thrust direction (=negative exhaust direction)
	double lsize;      ///<  exhaust length [m]
	double wsize;      ///<  exhaust width [m]
	double lofs;       ///<  longitudinal offset from engine [m]
	double modulate;   ///<  magnitude of random intensity variations (0..1)
	SURFHANDLE tex;    ///<  custom texture handle
	DWORD flags;       ///<  Bit flags (see \ref exhaustflag)
	UINT id;           ///<  reserved
} EXHAUSTSPEC;

/**
 * \brief Particle stream parameters
 * \note The following mapping methods (LEVELMAP) between stream
 *   level L and opacity a are supported:
 *   - LVL_FLAT: \f$ \alpha = \mathrm{const} \f$
 *   - LVL_LIN:  \f$ \alpha = L \f$
 *   - LVL_SQRT: \f$ \alpha = \sqrt(L) \f$
 *   - LVL_PLIN: \f$ \alpha = \left\lbrace \begin{array}{ll}
                     0 & \mathrm{if} L < L_\mathrm{min} \\
					 \frac{L-L_\mathrm{min}}{L_\mathrm{max}-L_\mathrm{min}} & \mathrm{if} L_\mathrm{min} \leq L \leq L_\mathrm{max} \\
					 1 & \mathrm{if} L > L_\mathrm{max}
					 \end{array} \right. \f$
 *   - LVL_PSQRT: \f$ \alpha = \left\lbrace \begin{array}{ll}
                     0 & \mathrm{if} L < L_\mathrm{min} \\
					 \sqrt{\frac{L-L_\mathrm{min}}{L_\mathrm{max}-L_\mathrm{min}}} & \mathrm{if} L_\mathrm{min} \leq L \leq L_\mathrm{max} \\
					 1 & \mathrm{if} L > L_\mathrm{max}
					 \end{array} \right. \f$
 */
typedef struct {
	DWORD flags;       ///<     streamspec bitflags
	double srcsize;    ///<     particle size at creation [m]
	double srcrate;    ///<     average particle creation rate [Hz]
	double v0;         ///<     emission velocity [m/s]
	double srcspread;  ///<     velocity spread during creation
	double lifetime;   ///<     average particle lifetime [s]
	double growthrate; ///<     particle growth rate [m/s]
	double atmslowdown;///<     slowdown rate in atmosphere
	/** \brief Particle lighting method */
	enum LTYPE {
		EMISSIVE,      ///<     emissive lighting (example: plasma stream)
		DIFFUSE        ///<     diffuse lighting (example: vapour stream)
	} ltype;           ///<     render lighting method
	/** \brief Mapping from level to alpha value (particle opacity) */
	enum LEVELMAP {
		LVL_FLAT,      ///<     constant (alpha independent of level)
		LVL_LIN,       ///<     linear mapping (alpha = level)
		LVL_SQRT,      ///<     square root mapping (alpha = sqrt(level)
		LVL_PLIN,      ///<     linear mapping in sub-range
		LVL_PSQRT      ///<     square-root mapping in sub-range
	} levelmap;        ///< mapping from level to alpha
	double lmin, lmax; ///<     min and max levels for level PLIN and PSQRT mapping types
	enum ATMSMAP { ATM_FLAT, ATM_PLIN, ATM_PLOG } atmsmap;    ///< mapping from atmospheric params to alpha
	double amin, amax; ///<     min and max densities for atms PLIN mapping
	SURFHANDLE tex;    ///<     particle texture handle (NULL for default)
} PARTICLESTREAMSPEC;

/** \brief Navigation transmitter data
 *
 * This structure contains both general data (transmitter type, channel,
 * output power and description string) and type-specific data.
 * To query type-specific data, first check the transmitter type, for example
 * \code
 * NAVDATA ndata;
 * oapiGetNavData (hNav, &ndata);
 * if (ndata.type == TRANSMITTER_ILS)
 *    approach_dir = ndata.ils.appdir;
 * \endcode
 * \note The power S<sub>0</sub> of a transmitter is defined in arbitrary units
 *   such that the signal S(r) = S<sub>0</sub>/r<sup>2</sup> drops to 1 at the
 *   maximum range r<sub>max</sub>, given a default receiver, i.e.
 *   S<sub>0</sub> = r<sup>2</sup><sub>max</sub>.
 * \sa oapiGetNavData
 */
typedef struct {
	// general data
	DWORD type;                ///< transmitter type id
	DWORD ch;                  ///< transmitter channel (0..639)
	double power;              ///< transmitter power [arbitrary units]
	const char *descr;         ///< pointer to transmitter description string
	// type-specific data
    union {
		struct {
			OBJHANDLE hPlanet; ///< associated planet
			double lng, lat;   ///< transmitter location [rad]
		} vor;
		struct {
			OBJHANDLE hBase;   ///< associated base
			int npad;          ///< pad number (>= 0)
		} vtol;
	    struct {
			OBJHANDLE hBase;   ///< associated base
		    double appdir;     ///< ILS approach direction [rad]
	    } ils;
		struct {
			OBJHANDLE hVessel; ///< associated vessel
			DOCKHANDLE hDock;  ///< associated docking port
		} ids;
		struct {
			OBJHANDLE hVessel; ///< associated vessel
		} xpdr;
	};
} NAVDATA;

/** \brief vessel beacon light parameters */
typedef struct {
	DWORD shape;       ///<   beacon shape identifier (see \ref beaconshape)
	VECTOR3 *pos;      ///<   pointer to position in vessel coordinates
	VECTOR3 *col;      ///<   pointer to beacon RGB colour
	double size;       ///<   beacon radius
	double falloff;    ///<   distance falloff parameter
	double period;     ///<   strobe period (0 for continuous)
	double duration;   ///<   strobe duration
	double tofs;       ///<   strobe time offset
	bool active;       ///<   beacon lit?
} BEACONLIGHTSPEC;

// ===========================================================================
/// \ingroup defines
/// \defgroup beaconshape Light beacon shape parameters
/// \sa BEACONLIGHTSPEC
// ===========================================================================
//@{
#define BEACONSHAPE_COMPACT 0 ///< compact beacon shape
#define BEACONSHAPE_DIFFUSE 1 ///< diffuse beacon shape
#define BEACONSHAPE_STAR    2 ///< star-shaped beacon
//@}

/**
 * \brief Vessel status parameters (version 1)
 *
 * Defines vessel status parameters at a given time. This is version 1 of
 * the vessel status interface. It is retained for backward compatibility,
 * but new modules should use VESSELSTATUS2 instead to exploit the latest
 * vessel capabilities such as individual thruster and propellant resource
 * settings.
 */
typedef struct {
	/// position relative to rbody in ecliptic frame [<b>m</b>]
	VECTOR3 rpos;

	/// velocity relative to rbody in ecliptic frame [<b>m/s</b>]
	VECTOR3 rvel;

	/// rotation velocity about principal axes in ecliptic frame [<b>rad/s</b>]
	VECTOR3 vrot;

	/// vessel orientation against ecliptic frame
	VECTOR3 arot;

	/// fuel level [0..1]
	double fuel;

	/// main/retro engine setting [-1..1]
	double eng_main;

	/// hover engine setting [0..1]
	double eng_hovr;

	/// handle of reference body
	OBJHANDLE rbody;

    /// handle of docking or landing target
	OBJHANDLE base;

    /// index of designated docking or landing port
	int port;

    /// \brief flight status indicator
	/// \note
	/// - 0=active (freeflight)
	/// - 1=inactive (landed)
	int status;

	/// \brief additional vector parameters
	/// \note
	/// - vdata[0]: contains landing paramters if status == 1:
	///   vdata[0].x = longitude, vdata[0].y = latitude, vdata[0].z = heading of landed vessel
	/// - vdata[1] - vdata[9]: not used
	VECTOR3 vdata[10];

	/// additional floating point parameters (not used)
	double  fdata[10];

	/// \brief additional integer and bitflag parameters
	///
	/// \par flag[0]&1:
	///   - 0: ingore eng_main and eng_hovr entries, do not change thruster settings
	///   - 1: set THGROUP_MAIN and THGROUP_RETRO thruster groups from eng_main, and THGROUP_HOVER from eng_hovr.
	/// \par flag[0]&2:
	///   - 0: ignore fuel level, do not change fuel levels
	///   - 1: set fuel level of first propellant resource from fuel
	/// \note flag[1] - flag[9]: not used
	DWORD   flag[10];
} VESSELSTATUS;


/**
 * \brief Vessel status parameters (version 2)
 * \details Defines vessel status parameters at a given time. This is version 2 of
 * the vessel status interface and replaces the earlier VESSELSTATUS
 * structure. Functions using VESSELSTATUS are still supported for backward
 * compatibility. \n
 * \note The version specification is an input parameter for all function calls
 *  (including GetStatus) and must be set by the user to tell Orbiter which interface to use.
 * \sa VESSEL::GetStatusEx
 */
typedef struct {
	/// interface version identifier (2)
	DWORD version;

	/**
	* \brief bit flags
	* \details The meaning of the bitflags in flag depends on whether the VESSELSTATUS2
	*  structure is used to get (GetStatus) or set (SetStatus) a vessel status. The
	*  following flags are currently defined:
	* \par flags:
	* - \c VS_FUELRESET
	*		- Get - not used
	*		- Set - reset all fuel levels to zero, independent of the fuel list.
	*		.
	* - \c VS_FUELLIST
	*		- Get - request a list of current fuel levels in fuel. The module is responsible
	*		  for deleting the list after use.
	*		- Set - set fuel levels for all resources listed in fuel.
	*		.
	* - \c VS_THRUSTRESET
	*		- Get - not used
	*		- Set - reset all thruster levels to zero, independent of the thruster list
	*		.
	* - \c VS_THRUSTLIST
	*		- Get - request a list of current thrust levels in thruster. The module is
	*		   responsible for deleting the list after use.
	*		- Set - set thrust levels for all thrusters listed in thruster.
	*		.
	* - \c VS_DOCKINFOLIST
	*		- Get - request a docking port status list in dockinfo. The module is
	*		  responsible for deleting the list after use.
	*		- Set - initialise docking status for all docking ports in dockinfo.
	* \sa VESSEL::GetStatusEx
	*/
	DWORD flag;

	/// handle of reference body
	OBJHANDLE rbody;

	/// handle of docking or landing target
	OBJHANDLE base;

	/// index of designated docking or landing port
	int port;

	/// \brief flight status indicator
	/// \note
	/// - 0=active (freeflight)
	/// - 1=inactive (landed)
	int status;

	/// position relative to reference body (rbody) in ecliptic frame [<b>m</b>]
	VECTOR3 rpos;

	/// velocity relative to reference body in ecliptic frame [<b>m/s</b>]
	VECTOR3 rvel;

	/// angular velocity around principal axes in ecliptic frame [<b>rad/s</b>]
	VECTOR3 vrot;

	/**
	* \brief vessel orientation against ecliptic frame
	* \details \b arot (\f$ \alpha, \beta, \gamma \f$) contains angles of rotation [rad]
	*  around \e x, \e y, \e z axes in ecliptic
	*  frame to produce this rotation matrix \b R for mapping from the vessel's local
	*  frame of reference to the global frame of reference:
	* \f[ R =
	*	\left[ \begin{array}{ccc} 1 & 0 & 0 \\	0 & \cos\alpha & \sin\alpha \\ 0 & -\sin\alpha & \cos\alpha \end{array} \right]
	*	\left[ \begin{array}{ccc} \cos\beta & 0 & -\sin\beta \\ 0 & 1 & 0 \\ \sin\beta & 0 & \cos\beta \end{array} \right]
	*	\left[ \begin{array}{ccc} \cos\gamma & \sin\gamma & 0 \\ -\sin\gamma & \cos\gamma & 0 \\ 0 & 0 & 1 \end{array} \right]
	* \f]
	* such that \b r<sub>global</sub> = \b R \b r<sub>local</sub> + \b p \n
	* where \b p is the vessel's global position.
	*/
	VECTOR3 arot;

	/**
	 * \brief longitude of vessel position in equatorial coordinates of rbody [rad]
	 * \note currently only defined if the vessel is landed (status=1)
	 */
	double surf_lng;

    /**
	 * \brief latitude of vessel position in equatorial coordinates of rbody [rad]
	 * \note currently only defined if the vessel is landed (status=1)
	 */
	double surf_lat;

	/**
	 * \brief vessel heading on the ground [rad]
	 * \note currently only defined if the vessel is landed (status=1)
	 */
	double surf_hdg;

	/// number of entries in the fuel list
	DWORD nfuel;

	/// propellant list
	struct FUELSPEC {
		DWORD idx;      ///< propellant index
		double level;   ///< propellant level
	} *fuel;

	/// number of entries in the thruster list
	DWORD nthruster;

	/// thruster definition list
	struct THRUSTSPEC {
		DWORD idx;      ///< thruster index
		double level;   ///< thruster level
	} *thruster;

	/// number of entries in the dockinfo list
	DWORD ndockinfo;

	/// dock info list
	struct DOCKINFOSPEC {
		DWORD idx;      ///< docking port index
		DWORD ridx;     ///< docking port index of docked vessel
		OBJHANDLE rvessel; ///< docked vessel
	} *dockinfo;

	/// transponder channel [0...640]
	DWORD xpdr;
} VESSELSTATUS2;

// VESSELSTATUSx bitflags
#define VS_FUELRESET    0x00000001 ///< set all propellant levels to zero
#define VS_FUELLIST     0x00000002 ///< list of propellant levels is provided
#define VS_THRUSTRESET  0x00000004 ///< set all thruster levels to zero
#define VS_THRUSTLIST   0x00000008 ///< list of thruster levels is provided
#define VS_DOCKINFOLIST 0x00000010 ///< list of docked objects is provided


/**
 * \brief Entry specification for selection list entry.
 */
typedef struct {
	char name[64];   ///< entry string
	DWORD flag;      ///< entry flags
} LISTENTRY;

/**
 * \brief Callback function for list entry selections.
 */
typedef bool (*Listentry_clbk)(char *name, DWORD idx, DWORD flag, void *usrdata);

/**
 * \ingroup defines
 * \defgroup listentryflag
 * \sa LISTENTRY
 */
//@
#define LISTENTRY_SUBITEM   0x01  ///< list entry has subitems
#define LISTENTRY_INACTIVE  0x02  ///< list entry can not be selected
#define LISTENTRY_SEPARATOR 0x04  ///< entry is followed by a separator
//@

#define LIST_UPENTRY 0x01         ///< list has parent list

/**
 * \ingroup defines
 * \defgroup listclbkflag
 * \sa LISTENTRY
 */
//@
#define LISTCLBK_CANCEL     0x00  ///< user cancelled the selection list
#define LISTCLBK_SELECT     0x01  ///< user selected an item
#define LISTCLBK_SUBITEM    0x02  ///< user steps down to subitem
#define LISTCLBK_UPLIST     0x03  ///< user steps up to parent list
//@

/**
 * \brief Context information for an Orbiter ingame help page.
 * \sa oapiOpenHelp
 */
typedef struct {
	char *helpfile;
	char *topic;
	char *toc;
	char *index;
} HELPCONTEXT;

typedef struct {
	char *name;
	void *parent;
	char *desc;
	void (*clbkFunc)(HINSTANCE,HWND);
} LP_EXTRAPRM;

#pragma pack(push,1)
/**
 * \brief This structure defines an affine mesh group transform
 *   (translation, rotation or scaling).
 * \sa VESSEL::MeshgroupTransform
 */
typedef struct {
	union {
		struct {
			VECTOR3 ref;   ///< rotation reference point
			VECTOR3 axis;  ///< rotation axis direction
			float angle;   ///< rotation angle [rad]
		} rotparam;
		struct {
			VECTOR3 shift; ///< translation vector
		} transparam;
		struct {
			VECTOR3 scale; ///< scaling factor
		} scaleparam;
	} P;
	int nmesh;             ///< mesh index (>= 0)
	int ngrp;              ///< group index (>= 0, or < 0 to indicate entire mesh)
	enum { TRANSLATE, ROTATE, SCALE } transform; ///< transformation flag
} MESHGROUP_TRANSFORM;
#pragma pack(pop)

// Animation component (obsolete)
typedef struct {
	UINT *grp;
	UINT ngrp;
	double state0;
	double state1;
	MESHGROUP_TRANSFORM trans;
} ANIMCOMP;


// Transformation base class
class MGROUP_TRANSFORM {
public:
	MGROUP_TRANSFORM () : mesh(0), grp(0), ngrp(0) {}
	MGROUP_TRANSFORM (UINT _mesh, UINT *_grp, UINT _ngrp)
		: mesh(_mesh), grp(_grp), ngrp(_ngrp) {}
	enum TYPE { NULLTRANSFORM, ROTATE, TRANSLATE, SCALE };
	virtual TYPE Type() const { return NULLTRANSFORM; }

	UINT mesh;
	UINT *grp;
	UINT ngrp;
};

// Rotation
class MGROUP_ROTATE: public MGROUP_TRANSFORM {
public:
	MGROUP_ROTATE (UINT _mesh, UINT *_grp, UINT _ngrp, const VECTOR3 &_ref, const VECTOR3 &_axis, float _angle)
		: MGROUP_TRANSFORM (_mesh, _grp, _ngrp), ref(_ref), axis(_axis), angle(_angle) {}
	TYPE Type() const { return ROTATE; }

	VECTOR3 ref;
	VECTOR3 axis;
	float angle;
};

// Translation
class MGROUP_TRANSLATE: public MGROUP_TRANSFORM {
public:
	MGROUP_TRANSLATE (UINT _mesh, UINT *_grp, UINT _ngrp, const VECTOR3 &_shift)
		: MGROUP_TRANSFORM (_mesh, _grp, _ngrp), shift(_shift) {}
	TYPE Type() const { return TRANSLATE; }

	VECTOR3 shift;
};

// Scaling
class MGROUP_SCALE: public MGROUP_TRANSFORM {
public:
	MGROUP_SCALE (UINT _mesh, UINT *_grp, UINT _ngrp, const VECTOR3 &_ref, const VECTOR3 &_scale)
		: MGROUP_TRANSFORM (_mesh, _grp, _ngrp), ref(_ref), scale(_scale) {}
	TYPE Type() const { return SCALE; }

	VECTOR3 ref;
	VECTOR3 scale;
};

/**
 * \brief Animation component definition
 *
 * Defines one component of an animation, including the mesh transformation,
 * the relative start and end points within the entire animation, and any
 * parent and child relationships with other animations.
 * \sa VESSEL::AddAnimationComponent
 */
struct ANIMATIONCOMP {
	double state0;			  ///< first end state
	double state1;			  ///< second end state
	MGROUP_TRANSFORM *trans;  ///< transformation
	ANIMATIONCOMP *parent;    ///< parent transformation
	ANIMATIONCOMP **children; ///< list of children
	UINT nchildren;           ///< number of children
};

/**
 * \brief Animation definition
 *
 * Defines a complete animation, including a list of components, the current
 * animation state, and the default state (as represented by the original mesh).
 */
struct ANIMATION {
	double defstate;          ///< default animation state in the mesh
	double state;             ///< current state
	UINT ncomp;               ///< number of components
	ANIMATIONCOMP **comp;     ///< list of components
};


/**
 * \ingroup defines
 * \defgroup animationflags Animation flags
 * \sa VESSEL::AddAnimationComponent
 */
//@{
#define LOCALVERTEXLIST ((UINT)(-1)) ///< flags animation component as explicit vertex list
#define MAKEGROUPARRAY(x) ((UINT*)x) ///< casts a vertex array into a group
//@}


typedef struct {
	RECT pos;
	int nbt_left, nbt_right;
	int bt_yofs, bt_ydist;
} MFDSPEC;

typedef struct {
	DWORD nmesh, ngroup;
} VCMFDSPEC;

#define MFD_SHOWMODELABELS 1

typedef struct {
	RECT pos;
	DWORD nmesh, ngroup;
	DWORD flag;
	int nbt1, nbt2;
	int bt_yofs, bt_ydist;
} EXTMFDSPEC;

typedef struct {
	DWORD nmesh, ngroup;
	VECTOR3 hudcnt;
	double size;
} VCHUDSPEC;

typedef struct {
	int W, H;
	int CX, CY;
	double Scale;
	int Markersize;
} HUDPAINTSPEC;

typedef struct {
	char *name;
	DWORD key;
	int (*msgproc)(UINT,UINT,WPARAM,LPARAM);
} MFDMODESPEC;

typedef struct {
	char *name;
	DWORD key;
	void *context;
	int (*msgproc)(UINT,UINT,WPARAM,LPARAM);
} MFDMODESPECEX;

typedef struct {
	int w, h;
	MFDMODESPECEX *spec;
} MFDMODEOPENSPEC;

#pragma pack(push,1)
typedef struct {
	const char *line1, *line2;
	char selchar;
} MFDBUTTONMENU;
#pragma pack(pop)


// ======================================================================
// Some helper functions
// ======================================================================

/**
 * \ingroup vec
 * \brief Vector composition
 *
 * Returns a vector composed of the three provided arguments
 * \param x x-component
 * \param y y-component
 * \param z z-component
 * \return vector defined as (x,y,z)
 */
inline VECTOR3 _V(double x, double y, double z)
{
	VECTOR3 vec = {x,y,z}; return vec;
}

/**
 * \ingroup vec
 * \brief Vector copy
 *
 * Copies the element values from the source to the target vector.
 * \param[out] a target vector
 * \param[in] b source vector
 */
inline void veccpy (VECTOR3 &a, const VECTOR3 &b)
{
	a.x = b.x;
	a.y = b.y;
	a.z = b.z;
}

/**
 * \ingroup vec
 * \brief Vector addition
 * \param a first vector operand
 * \param b second vector operand
 * \return Result of a+b.
 */
inline VECTOR3 operator+ (const VECTOR3 &a, const VECTOR3 &b)
{
	VECTOR3 c;
	c.x = a.x+b.x;
	c.y = a.y+b.y;
	c.z = a.z+b.z;
	return c;
}

/**
 * \ingroup vec
 * \brief Vector subtraction
 * \param a first vector operand
 * \param b second vector operand
 * \return Result of a-b.
 */
inline VECTOR3 operator- (const VECTOR3 &a, const VECTOR3 &b)
{
	VECTOR3 c;
	c.x = a.x-b.x;
	c.y = a.y-b.y;
	c.z = a.z-b.z;
	return c;
}

/**
 * \ingroup vec
 * \brief Multiplication of vector with scalar
 * \param a vector operand
 * \param f scalar operand
 * \return Result of element-wise a*f.
 */
inline VECTOR3 operator* (const VECTOR3 &a, const double f)
{
	VECTOR3 c;
	c.x = a.x*f;
	c.y = a.y*f;
	c.z = a.z*f;
	return c;
}

/**
 * \ingroup vec
 * \brief Division of vector by a scalar
 * \param a vector operand
 * \param f scalar operand
 * \return Result of element-wise a/f.
 */
inline VECTOR3 operator/ (const VECTOR3 &a, const double f)
{
	VECTOR3 c;
	c.x = a.x/f;
	c.y = a.y/f;
	c.z = a.z/f;
	return c;
}

/**
 * \ingroup vec
 * \brief Vector addition-assignment a += b
 * \param[in,out] a Left-hand vector operand
 * \param[in] b Right-hand vector operand
 * \return Replaces a with a+b and returns the result.
 */
inline VECTOR3 &operator+= (VECTOR3 &a, const VECTOR3 &b)
{
	a.x += b.x;
	a.y += b.y;
	a.z += b.z;
	return a;
}

/**
 * \ingroup vec
 * \brief Vector subtraction-assignment a -= b
 * \param[in,out] a Left-hand vector operand
 * \param[in] b Right-hand vector operand
 * \return Replaces a with a-b and returns the result.
 */
inline VECTOR3 &operator-= (VECTOR3 &a, const VECTOR3 &b)
{
	a.x -= b.x;
	a.y -= b.y;
	a.z -= b.z;
	return a;
}

/**
 * \ingroup vec
 * \brief Vector-scalar multiplication-assignment a *= f
 * \param[in,out] a Left-hand vector operand
 * \param[in] f Right hand scalar operand
 * \return Replaces a with element-wise a*f and returns the result.
 */
inline VECTOR3 &operator*= (VECTOR3 &a, const double f)
{
	a.x *= f;
	a.y *= f;
	a.z *= f;
	return a;
}

/**
 * \ingroup vec
 * \brief Vector-scalar division-assignment a /= f
 * \param[in,out] a Left-hand vector operand
 * \param[in] f Right-hand scalar operand
 * \return Replaces a with element-wise a/f and returns the result.
 */
inline VECTOR3 &operator/= (VECTOR3 &a, const double f)
{
	a.x /= f;
	a.y /= f;
	a.z /= f;
	return a;
}

/**
 * \ingroup vec
 * \brief Vector unary minus -a
 * \param[in] a Vector operand
 * \return Negative vector (-a.x, -a.y, -a.z)
 */
inline VECTOR3 operator- (const VECTOR3 &a)
{
	VECTOR3 c;
	c.x = -a.x;
	c.y = -a.y;
	c.z = -a.z;
	return c;
}

/**
 * \ingroup vec
 * \brief Scalar (inner, dot) product of two vectors
 * \param[in] a First vector operand
 * \param[in] b Second vector operand
 * \return Scalar product <b>ab</b>
 */
inline double dotp (const VECTOR3 &a, const VECTOR3 &b)
{
	return a.x*b.x + a.y*b.y + a.z*b.z;
}

/**
 * \ingroup vec
 * \brief Vector (cross) product of two vectors
 * \param[in] a First vector operand
 * \param[in] b Second vector operand
 * \return Vector product <b>a</b>x<b>b</b>
 */
inline VECTOR3 crossp (const VECTOR3 &a, const VECTOR3 &b)
{
	return _V(a.y*b.z - b.y*a.z, a.z*b.x - b.z*a.x, a.x*b.y - b.x*a.y);
}

/**
 * \ingroup vec
 * \brief Length (L2-norm) of a vector
 * \param a Vector operand
 * \return Vector norm |<b>a</b>|<sub>2</sub>
 */
inline double length (const VECTOR3 &a)
{
	return sqrt (a.x*a.x + a.y*a.y + a.z*a.z);
}

/**
 * \ingroup vec
 * \brief Distance between two points
 * \param[in] a First point
 * \param[in] b Second point
 * \return Distance between a and b
 */
inline double dist (const VECTOR3 &a, const VECTOR3 &b)
{
	return length (a-b);
}

/**
 * \ingroup vec
 * \brief Normalise a vector
 *
 * Resizes the argument vector to length 1.
 * \param[in,out] a Vector argument
 * \note The length of a must be greater than 0.
 */
inline void normalise (VECTOR3 &a)
{
	a /= length(a);
}

/**
 * \ingroup vec
 * \brief Returns normalised vector
 *
 * Returns a vector of length 1 with the same direction
 * as the argument vector.
 * \param[in] a Vector argument
 * \return Normalised vector.
 * \note The length of a must be greater than 0.
 */
inline VECTOR3 unit (const VECTOR3 &a)
{
	return a / length(a);
}

/**
 * \ingroup vec
 * \brief Matrix composition
 *
 * Returns a matrix composed of the provided elements.
 * \return
 * \f$
 *  \left(\begin{array}{ccc}
 *  m_{11} & m_{12} & m_{13} \\
 *  m_{21} & m_{22} & m_{23} \\
 *  m_{31} & m_{32} & m_{33}
 *  \end{array}\right)
 * \f$
 */
inline MATRIX3 _M(double m11, double m12, double m13,
				  double m21, double m22, double m23,
				  double m31, double m32, double m33)
{
	MATRIX3 mat = {m11,m12,m13,  m21,m22,m23,  m31,m32,m33};
	return mat;
}

/**
 * \ingroup vec
 * \brief Returns the identity matrix
 */
inline MATRIX3 identity ()
{
	static MATRIX3 mat = {1,0,0, 0,1,0, 0,0,1};
	return mat;
}

/**
 * \ingroup vec
 * \brief Outer product of two vectors
 * \param[in] a First vector operand
 * \param[in] b Second vector operand
 * \return Outer product <b>a</b><b>b</b><sup>T</sup>, where
 * <b>a</b> and <b>b</b> represent column vectors.
 */
inline MATRIX3 outerp (const VECTOR3 &a, const VECTOR3 &b)
{
	return _M(a.x*b.x, a.x*b.y, a.x*b.z,
		      a.y*b.x, a.y*b.y, a.y*b.z,
			  a.z*b.x, a.z*b.y, a.z*b.z);
}

/**
 * \ingroup vec
 * \brief Sum of matrix and scalar.
 * \param[in] A Matrix operand (left)
 * \param[in] s scalar operand (right)
 * \return A+s (element-wise sum of A and s)
 */
inline MATRIX3 operator+ (const MATRIX3 &A, double s)
{
	MATRIX3 mat = {A.m11+s, A.m12+s, A.m13+s,
		           A.m21+s, A.m22+s, A.m23+s,
				   A.m31+s, A.m32+s, A.m33+s};
	return mat;
}

/**
 * \ingroup vec
 * \brief Difference of matrix and scalar.
 * \param[in] A Matrix operand (left)
 * \param[in] s scalar operand (right)
 * \return A-s (element-wise difference of A and s)
 */
inline MATRIX3 operator- (const MATRIX3 &A, double s)
{
	MATRIX3 mat = {A.m11-s, A.m12-s, A.m13-s,
		           A.m21-s, A.m22-s, A.m23-s,
				   A.m31-s, A.m32-s, A.m33-s};
	return mat;
}

/**
 * \ingroup vec
 * \brief Product of matrix and scalar.
 * \param[in] A Matrix operand (left)
 * \param[in] s scalar operand (right)
 * \return A*s (element-wise product of A and s)
 */
inline MATRIX3 operator* (const MATRIX3 &A, double s)
{
	MATRIX3 mat = {A.m11*s, A.m12*s, A.m13*s,
		           A.m21*s, A.m22*s, A.m23*s,
				   A.m31*s, A.m32*s, A.m33*s};
	return mat;
}

/**
 * \ingroup vec
 * \brief Quotient of matrix and scalar.
 * \param[in] A Matrix operand (left)
 * \param[in] s scalar operand (right)
 * \return A/s (element-wise quotient of A and s)
 * \note s != 0 is required.
 */
inline MATRIX3 operator/ (const MATRIX3 &A, double s)
{
	MATRIX3 mat = {A.m11/s, A.m12/s, A.m13/s,
		           A.m21/s, A.m22/s, A.m23/s,
				   A.m31/s, A.m32/s, A.m33/s};
	return mat;
}

/**
 * \ingroup vec
 * \brief Matrix-scalar product-assignment A *= s
 * \param[in] A Matrix operand (left)
 * \param[in] s scalar operand (right)
 * \return Replaces A with element-wise product A*s and returns the result.
 */
inline MATRIX3 &operator*= (MATRIX3 &A, double s)
{
	for (int i = 0; i < 9; i++) A.data[i] *= s;
	return A;
}

/**
 * \ingroup vec
 * \brief Matrix-scalar division-assignment A /= s
 * \param[in] A Matrix operand (left)
 * \param[in] s scalar operand (right)
 * \return Replaces A with element-wise quotient A/s and returns the result.
 * \note s != 0 is required.
 */
inline MATRIX3 &operator/= (MATRIX3 &A, double s)
{
	for (int i = 0; i < 9; i++) A.data[i] /= s;
	return A;
}

/**
 * \ingroup vec
 * \brief Matrix-vector multiplication
 * \param[in] A matrix operand
 * \param[in] b vector operand
 * \return Result of <b>Ab</b>
 */
inline VECTOR3 mul (const MATRIX3 &A, const VECTOR3 &b)
{
	return _V (
		A.m11*b.x + A.m12*b.y + A.m13*b.z,
		A.m21*b.x + A.m22*b.y + A.m23*b.z,
		A.m31*b.x + A.m32*b.y + A.m33*b.z);
}

/**
 * \ingroup vec
 * \brief Matrix transpose-vector multiplication
 * \param[in] A matrix operand
 * \param[in] b vector operand
 * \return Result of <b>A</b><sup>T</sup><b>b</b>
 */
inline VECTOR3 tmul (const MATRIX3 &A, const VECTOR3 &b)
{
	return _V (
		A.m11*b.x + A.m21*b.y + A.m31*b.z,
		A.m12*b.x + A.m22*b.y + A.m32*b.z,
		A.m13*b.x + A.m23*b.y + A.m33*b.z);
}

/**
 * \ingroup vec
 * \brief Matrix-matrix multiplication
 * \param[in] A First matrix operand
 * \param[in] B Second matrix operand
 * \return Result of <b>AB</b>
 */
inline MATRIX3 mul (const MATRIX3 &A, const MATRIX3 &B)
{
	MATRIX3 mat = {
		A.m11*B.m11 + A.m12*B.m21 + A.m13*B.m31, A.m11*B.m12 + A.m12*B.m22 + A.m13*B.m32, A.m11*B.m13 + A.m12*B.m23 + A.m13*B.m33,
		A.m21*B.m11 + A.m22*B.m21 + A.m23*B.m31, A.m21*B.m12 + A.m22*B.m22 + A.m23*B.m32, A.m21*B.m13 + A.m22*B.m23 + A.m23*B.m33,
		A.m31*B.m11 + A.m32*B.m21 + A.m33*B.m31, A.m31*B.m12 + A.m32*B.m22 + A.m33*B.m32, A.m31*B.m13 + A.m32*B.m23 + A.m33*B.m33
	};
	return mat;
}

inline RECT _R (int left, int top, int right, int bottom)
{
	RECT r = { left, top, right, bottom }; return r;
}

inline VECTOR3 POINTERTOREF (VECTOR3 *p)
{
	VECTOR3 v;
	v.x = std::numeric_limits<double>::max();            // flag
	*((VECTOR3**)&v.z) = p;   // address
	v.z = 0.0;
	return v;
}


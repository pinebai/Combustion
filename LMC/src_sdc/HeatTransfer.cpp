//
// "Divu_Type" means S, where divergence U = S
// "Dsdt_Type" means pd S/pd t, where S is as above
// "RhoYchemProd_Type" means -omega_l/rho, i.e., the mass rate of decrease of species l due
//             to kinetics divided by rho
//
#include <winstd.H>

#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cfloat>
#include <fstream>
#include <vector>

#include <Geometry.H>
#include <BoxDomain.H>
#include <ParmParse.H>
#include <ErrorList.H>
#include <HeatTransfer.H>
#include <HEATTRANSFER_F.H>
#include <DIFFUSION_F.H>
#include <MultiGrid.H>
#include <ArrayLim.H>
#include <SPACE.H>
#include <Interpolater.H>
#include <ccse-mpi.H>
#include <Utility.H>

#if defined(BL_USE_NEWMECH) || defined(BL_USE_VELOCITY)
#include <DataServices.H>
#include <AmrData.H>
#endif

#include <PROB_F.H>
#include <NAVIERSTOKES_F.H>
#include <VISCOPERATOR_F.H>
#include <DERIVE_F.H>


static Box stripBox; // used for debugging


using std::cout;
using std::endl;
using std::cerr;

#define DEF_LIMITS(fab,fabdat,fablo,fabhi)   \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITS(fab,fabdat,fablo,fabhi)  \
const int* fablo = (fab).loVect();           \
const int* fabhi = (fab).hiVect();           \
const Real* fabdat = (fab).dataPtr();

#define DEF_CLIMITSCOMP(fab,fabdat,fablo,fabhi,comp)  \
const int* fablo = (fab).loVect();                    \
const int* fabhi = (fab).hiVect();                    \
const Real* fabdat = (fab).dataPtr(comp);

#ifdef BL_USE_FLOAT
#  define Real_MIN FLT_MIN
#  define Real_MAX FLT_MAX
#else
#  define Real_MIN DBL_MIN
#  define Real_MAX DBL_MAX
#endif

#define GEOM_GROW   1
#define HYP_GROW    3
#define PRESS_GROW  1
#define DIVU_GROW   1
#define DSDT_GROW   1
#define bogus_value 1.e20
#define DQRAD_GROW  1
#define YDOT_GROW   1
#define HYPF_GROW   1

const int LinOp_grow = 1;

static const std::string typical_values_filename("typical_values.fab");

namespace
{
    bool initialized = false;
}
//
// Set all default values in Initialize()!!!
//
namespace
{
    std::set<std::string> ShowMF_Sets;
    std::string           ShowMF_Dir;
    bool                  ShowMF_Verbose;
    bool                  ShowMF_Check_Nans;
    FABio::Format         ShowMF_Fab_Format;
    bool                  do_not_use_funccount;
    bool                  do_active_control;
    Real                  crse_dt;
}

int  HeatTransfer::num_divu_iters;
int  HeatTransfer::init_once_done;
int  HeatTransfer::do_OT_radiation;
int  HeatTransfer::do_heat_sink;
int  HeatTransfer::RhoH;
int  HeatTransfer::do_diffuse_sync;
int  HeatTransfer::dpdt_option;
int  HeatTransfer::RhoYdot_Type;
int  HeatTransfer::FuncCount_Type;
int  HeatTransfer::divu_ceiling;
Real HeatTransfer::divu_dt_factor;
Real HeatTransfer::min_rho_divu_ceiling;
int  HeatTransfer::have_trac;
int  HeatTransfer::have_rhort;
int  HeatTransfer::Trac;
int  HeatTransfer::RhoRT;
int  HeatTransfer::first_spec;
int  HeatTransfer::last_spec;
int  HeatTransfer::nspecies;
int  HeatTransfer::floor_species;
int  HeatTransfer::do_set_rho_to_species_sum;
Real HeatTransfer::rgas;
Real HeatTransfer::prandtl;
Real HeatTransfer::schmidt;
Real HeatTransfer::constant_mu_val;
Real HeatTransfer::constant_rhoD_val;
Real HeatTransfer::constant_lambda_val;
int  HeatTransfer::unity_Le;
Real HeatTransfer::htt_tempmin;
Real HeatTransfer::htt_tempmax;
Real HeatTransfer::htt_hmixTYP;
int  HeatTransfer::zeroBndryVisc;
int  HeatTransfer::do_add_nonunityLe_corr_to_rhoh_adv_flux;
int  HeatTransfer::do_check_divudt;
int  HeatTransfer::hack_nochem;
int  HeatTransfer::hack_nospecdiff;
int  HeatTransfer::hack_nomcddsync;
int  HeatTransfer::hack_noavgdivu;
Real HeatTransfer::trac_diff_coef;
Real HeatTransfer::P1atm_MKS;
bool HeatTransfer::plot_reactions;
bool HeatTransfer::plot_consumption;
bool HeatTransfer::plot_heat_release;
int  HeatTransfer::do_mcdd;
int  HeatTransfer::mcdd_NitersMAX;
Real HeatTransfer::mcdd_relax_factor_T;
Real HeatTransfer::mcdd_relax_factor_Y;
int  HeatTransfer::mcdd_mgLevelsMAX;
int  HeatTransfer::mcdd_nub;
int  HeatTransfer::mcdd_numcycles;
int  HeatTransfer::mcdd_verbose;
Real HeatTransfer::mcdd_res_nu1_rtol;
Real HeatTransfer::mcdd_res_nu2_rtol;
Real HeatTransfer::mcdd_res_redux_tol;
Real HeatTransfer::mcdd_res_abs_tol;
Real HeatTransfer::mcdd_stalled_tol;
Real HeatTransfer::mcdd_advance_temp;
Real HeatTransfer::new_T_threshold;
int  HeatTransfer::nGrowAdvForcing=1;

std::string                                HeatTransfer::turbFile;
ChemDriver*                                HeatTransfer::chemSolve;
std::map<std::string, Array<std::string> > HeatTransfer::auxDiag_names;
std::string                                HeatTransfer::mcdd_transport_model;

Array<int>  HeatTransfer::mcdd_nu1;
Array<int>  HeatTransfer::mcdd_nu2;
Array<Real> HeatTransfer::typical_values;

// Turns out to be handy for debugging
#if 1
static
void dump(const FArrayBox& fab, int comp=-1)
{

#if BL_SPACEDIM == 2
    int sComp = comp < 0 ? 0 : comp;
    int nComp = comp < 0 ? fab.nComp() : 1;

    for (int j=62; j<=62; ++j) {
        IntVect iv(j,18);
        cout << j << " ";
        for (int n=sComp; n<sComp+nComp; ++n) {
            cout << fab(iv,n) << " ";
        }
        cout << endl;
    }
#endif

}
static
void dump(const MultiFab& mf, int comp=-1)
{
    dump(mf[0],comp);
}
#endif

///////////////////////////////
// SDC Stuff

// these can be set in the inputs file
int HeatTransfer::sdc_iterMAX;

#ifdef PARTICLES

namespace
{
    //
    // Name of subdirectory in chk???? holding checkpointed particles.
    // 
    const std::string the_ht_particle_file_name("Particles");
    //
    // There's really only one of these.
    //
    HTParticleContainer* HTPC = 0;

    std::string      timestamp_dir;
    std::vector<int> timestamp_indices;
    std::string      particle_init_file;
    std::string      particle_restart_file;
    std::string      particle_output_file;
    bool             restart_from_nonparticle_chkfile;
    int              pverbose;
    //
    // We want to call this routine on exit to clean up particles.
    //
    void RemoveParticles ()
    {
        delete HTPC;
        HTPC = 0;
    }
}
//
// In case someone outside of HeatTransfer needs a handle on the particles.
//
HTParticleContainer* HeatTransfer::theHTPC () { return HTPC; }

#endif /*PARTICLES*/

static
std::string to_upper(const std::string& s)
{
  std::string rtn = s;
  for (unsigned int i=0; i<rtn.length(); i++)
  {
      rtn[i] = toupper(rtn[i]);
  }
  return rtn;
}

static std::map<std::string,FABio::Format> ShowMF_Fab_Format_map;


void
HeatTransfer::Initialize ()
{
    if (initialized) return;

    NavierStokes::Initialize();
    //
    // Set all default values here!!!
    //
    ShowMF_Fab_Format_map["ASCII"] = FABio::FAB_ASCII;
    ShowMF_Fab_Format_map["IEEE"] = FABio::FAB_IEEE;
    ShowMF_Fab_Format_map["NATIVE"] = FABio::FAB_NATIVE;
    ShowMF_Fab_Format_map["8BIT"] = FABio::FAB_8BIT;
    ShowMF_Fab_Format_map["IEEE_32"] = FABio::FAB_IEEE_32;
    ShowMF_Verbose       = true;
    ShowMF_Check_Nans    = true;
    ShowMF_Fab_Format    = ShowMF_Fab_Format_map["ASCII"];
    do_not_use_funccount = false;
    do_active_control    = false;
    crse_dt              = -1;
    
    HeatTransfer::num_divu_iters            = 1;
    HeatTransfer::init_once_done            = 0;
    HeatTransfer::do_OT_radiation           = 0;
    HeatTransfer::do_heat_sink              = 0;
    HeatTransfer::RhoH                      = -1;
    HeatTransfer::do_diffuse_sync           = 1;
    HeatTransfer::dpdt_option               = 2;
    HeatTransfer::RhoYdot_Type                 = -1;
    HeatTransfer::FuncCount_Type            = -1;
    HeatTransfer::divu_ceiling              = 0;
    HeatTransfer::divu_dt_factor            = .5;
    HeatTransfer::min_rho_divu_ceiling      = -1.e20;
    HeatTransfer::have_trac                 = 0;
    HeatTransfer::have_rhort                = 0;
    HeatTransfer::Trac                      = -1;
    HeatTransfer::RhoRT                     = -1;
    HeatTransfer::first_spec                = -1;
    HeatTransfer::last_spec                 = -2;
    HeatTransfer::nspecies                  = 0;
    HeatTransfer::floor_species             = 0;
    HeatTransfer::do_set_rho_to_species_sum = 1;
    HeatTransfer::rgas                      = -1.0;
    HeatTransfer::prandtl                   = .7;
    HeatTransfer::schmidt                   = .7;
    HeatTransfer::constant_mu_val           = -1;
    HeatTransfer::constant_rhoD_val         = -1;
    HeatTransfer::constant_lambda_val       = -1;
    HeatTransfer::unity_Le                  = 1;
    HeatTransfer::htt_tempmin               = 298.0;
    HeatTransfer::htt_tempmax               = 40000.;
    HeatTransfer::htt_hmixTYP               = -1.;
    HeatTransfer::zeroBndryVisc             = 0;
    HeatTransfer::chemSolve                 = 0;
    HeatTransfer::do_check_divudt           = 1;
    HeatTransfer::hack_nochem               = 0;
    HeatTransfer::hack_nospecdiff           = 0;
    HeatTransfer::hack_nomcddsync           = 1;
    HeatTransfer::hack_noavgdivu            = 0;
    HeatTransfer::trac_diff_coef            = 0.0;
    HeatTransfer::P1atm_MKS                 = -1.0;
    HeatTransfer::turbFile                  = "";
    HeatTransfer::plot_reactions            = false;
    HeatTransfer::plot_consumption          = true;
    HeatTransfer::plot_heat_release         = true;
    HeatTransfer::do_mcdd                   = 0;
    HeatTransfer::mcdd_transport_model      = "";
    HeatTransfer::mcdd_NitersMAX            = 100;
    HeatTransfer::mcdd_relax_factor_T       = 1.0;
    HeatTransfer::mcdd_relax_factor_Y       = 1.0;
    HeatTransfer::mcdd_mgLevelsMAX          = -1;
    HeatTransfer::mcdd_nub                  = -1;
    HeatTransfer::mcdd_numcycles            = 10;
    HeatTransfer::mcdd_verbose              = 1;
    HeatTransfer::mcdd_res_nu1_rtol         = 1.e-8;
    HeatTransfer::mcdd_res_nu2_rtol         = 1.e-8;
    HeatTransfer::mcdd_res_redux_tol        = 1.e-8;
    HeatTransfer::mcdd_res_abs_tol          = 1.e-8;
    HeatTransfer::mcdd_stalled_tol          = 1.e-8;
    HeatTransfer::mcdd_advance_temp         = 1;
    HeatTransfer::new_T_threshold           = -1;  // On new AMR level, max change in lower bound for T, not used if <=0

    HeatTransfer::do_add_nonunityLe_corr_to_rhoh_adv_flux = 1;

    HeatTransfer::sdc_iterMAX               = 1;

#ifdef PARTICLES
    timestamp_dir                    = "Timestamps";
    restart_from_nonparticle_chkfile = false;
    pverbose                         = 2;
#endif /*PARTICLES*/

    ParmParse pp("ns");

    pp.query("do_diffuse_sync",do_diffuse_sync);
    BL_ASSERT(do_diffuse_sync == 0 || do_diffuse_sync == 1);
    pp.query("dpdt_option",dpdt_option);
    BL_ASSERT(dpdt_option >= 0 && dpdt_option <= 2);
    pp.query("do_active_control",do_active_control);

    verbose = pp.contains("v");

    pp.query("divu_ceiling",divu_ceiling);
    BL_ASSERT(divu_ceiling >= 0 && divu_ceiling <= 3);
    pp.query("divu_dt_factor",divu_dt_factor);
    BL_ASSERT(divu_dt_factor>0 && divu_dt_factor <= 1.0);
    pp.query("min_rho_divu_ceiling",min_rho_divu_ceiling);
    if (divu_ceiling) BL_ASSERT(min_rho_divu_ceiling >= 0.0);

    pp.query("htt_tempmin",htt_tempmin);
    pp.query("htt_tempmax",htt_tempmax);

    pp.query("floor_species",floor_species);
    BL_ASSERT(floor_species == 0 || floor_species == 1);

    pp.query("do_set_rho_to_species_sum",do_set_rho_to_species_sum);

    pp.query("num_divu_iters",num_divu_iters);

    pp.query("do_not_use_funccount",do_not_use_funccount);

    pp.query("schmidt",schmidt);
    pp.query("prandtl",prandtl);
    pp.query("unity_Le",unity_Le);
    unity_Le = unity_Le ? 1 : 0;
    if (unity_Le)
    {
	schmidt = prandtl;
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: Le=1, setting Sc = Pr" << '\n';
    }

    pp.query("sdc_iterMAX",sdc_iterMAX);

    pp.query("constant_mu_val",constant_mu_val);
    pp.query("constant_rhoD_val",constant_rhoD_val);
    pp.query("constant_lambda_val",constant_lambda_val);
    if (constant_mu_val != -1)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: using constant_mu_val = " 
                      << constant_mu_val << '\n';
    }
    if (constant_rhoD_val != -1)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: using constant_rhoD_val = " 
                      << constant_rhoD_val << '\n';
    }
    if (constant_lambda_val != -1)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::read_params: using constant_lambda_val = " 
                      << constant_lambda_val << '\n';
    }

    pp.query("do_add_nonunityLe_corr_to_rhoh_adv_flux", do_add_nonunityLe_corr_to_rhoh_adv_flux);
    pp.query("hack_nochem",hack_nochem);
    pp.query("hack_nospecdiff",hack_nospecdiff);
    pp.query("hack_nomcddsync",hack_nomcddsync);
    pp.query("hack_noavgdivu",hack_noavgdivu);
    pp.query("do_check_divudt",do_check_divudt);

    pp.query("do_OT_radiation",do_OT_radiation);
    do_OT_radiation = (do_OT_radiation ? 1 : 0);
    pp.query("do_heat_sink",do_heat_sink);
    do_heat_sink = (do_heat_sink ? 1 : 0);

    chemSolve = new ChemDriver();

    pp.query("turbFile",turbFile);

    pp.query("zeroBndryVisc",zeroBndryVisc);
    //
    // Set variability/visc for velocities.
    //
    if (variable_vel_visc != 1)
      BoxLib::Error("HeatTransfer::read_params() -- must use variable viscosity");
    //
    // Set variability/visc for tracer
    //
    if (variable_scal_diff != 1)
      BoxLib::Error("HeatTransfer::read_params() -- must use variable scalar diffusivity");
    //
    // Read in scalar value and use it as tracer.
    //
    pp.query("scal_diff_coefs",trac_diff_coef);

    for (int i = 0; i < visc_coef.size(); i++)
        visc_coef[i] = bogus_value;

    // Useful for debugging
    ParmParse pproot;
    if (int nsv=pproot.countval("ShowMF_Sets"))
    {
        Array<std::string> ShowMF_set_names(nsv);
        pproot.getarr("ShowMF_Sets",ShowMF_set_names);
        for (int i=0; i<nsv; ++i) {
            ShowMF_Sets.insert(ShowMF_set_names[i]);
        }
        ShowMF_Dir="."; pproot.query("ShowMF_Dir",ShowMF_Dir);
        pproot.query("ShowMF_Verbose",ShowMF_Verbose);
        pproot.query("ShowMF_Check_Nans",ShowMF_Check_Nans);
        std::string format="NATIVE"; pproot.query("ShowMF_Fab_Format",format);
        if (ShowMF_Fab_Format_map.count(to_upper(format)) == 0) {
            BoxLib::Abort("Unknown FABio format label");
        }
        ShowMF_Fab_Format = ShowMF_Fab_Format_map[format];

        if (ShowMF_Verbose>0 && ShowMF_set_names.size()>0 && ParallelDescriptor::IOProcessor()) {
            cout << "   ******************************  Debug: ShowMF_Sets: ";
            for (int i=0; i<ShowMF_set_names.size(); ++i) {
                cout << ShowMF_set_names[i] << " ";
            }
            cout << endl;
        }

    }

    pp.query("do_mcdd",do_mcdd);
    if (do_mcdd) {
        DDOp::set_chem_driver(getChemSolve());
        pp.get("mcdd_transport_model",mcdd_transport_model);
        if (mcdd_transport_model=="full") {
            DDOp::set_transport_model(DDOp::DD_Model_Full);
        }
        else if (mcdd_transport_model=="mix") {
            DDOp::set_transport_model(DDOp::DD_Model_MixAvg);
        }
        else
        {
            BoxLib::Abort("valid models: full, mix");
        }
        pp.query("mcdd_NitersMAX",mcdd_NitersMAX);
        pp.query("mcdd_relax_factor_T",mcdd_relax_factor_T);
        pp.query("mcdd_relax_factor_Y",mcdd_relax_factor_Y);
        pp.query("mcdd_numcycles",mcdd_numcycles);
        if (mcdd_numcycles<=0) mcdd_numcycles=1; // 0 is not valid, assume user wants one cycle
        pp.query("mcdd_mgLevelsMAX",mcdd_mgLevelsMAX);
        if (mcdd_mgLevelsMAX==0) mcdd_mgLevelsMAX=1; // 0 is not valid, assume user wants no additional levels
        DDOp::set_mgLevelsMAX(mcdd_mgLevelsMAX);

        int npre = pp.countval("mcdd_nu1");
        if (npre>0) {
            mcdd_nu1.resize(npre);
            pp.getarr("mcdd_nu1",mcdd_nu1,0,npre);
        }
        int max_nvals = (mcdd_mgLevelsMAX<0  ?  100  :  mcdd_mgLevelsMAX);
        mcdd_nu1.resize(max_nvals);
        for (int i=npre; i<max_nvals; ++i) {
            mcdd_nu1[i] = mcdd_nu1[npre-1];
        }
        int npost = pp.countval("mcdd_nu2");
        if (npost>0) {
            mcdd_nu2.resize(npost);
            pp.getarr("mcdd_nu2",mcdd_nu2,0,npost);
        }
        mcdd_nu2.resize(max_nvals);
        for (int i=npost; i<max_nvals; ++i) {
            mcdd_nu2[i] = mcdd_nu2[npost-1];
        }
        pp.query("mcdd_nub",mcdd_nub);
        pp.query("mcdd_res_nu1_rtol",mcdd_res_nu1_rtol);
        pp.query("mcdd_res_nu2_rtol",mcdd_res_nu2_rtol);
        pp.query("mcdd_res_redux_tol",mcdd_res_redux_tol);
        pp.query("mcdd_res_abs_tol",mcdd_res_abs_tol);
        pp.query("mcdd_stalled_tol",mcdd_stalled_tol);
        pp.query("mcdd_advance_temp",mcdd_advance_temp);

        pp.query("mcdd_verbose",mcdd_verbose);
    }

#ifdef PARTICLES
    //
    // Some particle stuff.
    //
    ParmParse ppp("particles");
    //
    // The directory in which to store timestamp files.
    //
    ppp.query("timestamp_dir", timestamp_dir);
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(timestamp_dir, 0755))
            BoxLib::CreateDirectoryFailed(timestamp_dir);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (int nc = ppp.countval("timestamp_indices"))
    {
        timestamp_indices.resize(nc);

        ppp.getarr("timestamp_indices", timestamp_indices, 0, nc);
    }

    ppp.query("pverbose",pverbose);
    //
    // Used in initData() on startup to read in a file of particles.
    //
    ppp.query("particle_init_file", particle_init_file);
    //
    // Used in post_restart() to read in a file of particles.
    //
    ppp.query("particle_restart_file", particle_restart_file);
    //
    // This must be true the first time you try to restart from a checkpoint
    // that was written with USE_PARTICLES=FALSE; i.e. one that doesn't have
    // the particle checkpoint stuff (even if there are no active particles).
    // Otherwise the code will fail when trying to read the checkpointed particles.
    //
    ppp.query("restart_from_nonparticle_chkfile", restart_from_nonparticle_chkfile);
    //
    // Used in post_restart() to write out the file of particles.
    //
    ppp.query("particle_output_file", particle_output_file);
#endif
        
    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nDumping ParmParse table:\n \n";
        ParmParse::dumpTable(std::cout);
        std::cout << "\n... done dumping ParmParse table.\n" << '\n';
    }

    BoxLib::ExecOnFinalize(HeatTransfer::Finalize);

    initialized = true;
}

void
HeatTransfer::Finalize ()
{
    initialized = false;
}

static
void
FabMinMax (FArrayBox& fab,
           const Box& box,
           Real       fmin,
           Real       fmax,
           int        sComp,
           int        nComp)
{
    BL_ASSERT(fab.box().contains(box));
    BL_ASSERT(sComp + nComp <= fab.nComp());

    const int* lo     = box.loVect();
    const int* hi     = box.hiVect();
    Real*      fabdat = fab.dataPtr(sComp);
    
    FORT_FABMINMAX(lo, hi,
                   fabdat, ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
                   &fmin, &fmax, &nComp);
}

static
std::ostream&
levWrite(std::ostream& os, int level, std::string& message)
{
    if (ParallelDescriptor::IOProcessor())
    {
        for (int lev=0; lev<=level; ++lev) {
            //os << '\t';
            os << "  ";
        }
        os << message;
    }
    return os;
}

Box
getStrip(const Geometry& geom)
{
    const Box& box = geom.Domain();
    IntVect be = box.bigEnd();
    IntVect se = box.smallEnd();
    se[0] = (int) 0.5*(se[0]+be[0]);
    be[0] = se[0];
    return Box(se,be);
}

void
showMFsub(const std::string&   mySet,
          const MultiFab&      mf,
          const Box&           box,
          const std::string&   name,
          int                  lev = -1,
          int                  iter = -1) // Default value = no append 2nd integer
{
    if (ShowMF_Sets.count(mySet)>0)
    {
        const FABio::Format saved_format = FArrayBox::getFormat();
        FArrayBox::setFormat(ShowMF_Fab_Format);
        std::string DebugDir(ShowMF_Dir);
        if (ParallelDescriptor::IOProcessor())
            if (!BoxLib::UtilCreateDirectory(DebugDir, 0755))
                BoxLib::CreateDirectoryFailed(DebugDir);
        ParallelDescriptor::Barrier();

        std::string junkname = name;
        if (lev>=0) {
            junkname = BoxLib::Concatenate(junkname+"_",lev,1);
        }
        if (iter>=0) {
            junkname = BoxLib::Concatenate(junkname+"_",iter,1);
        }
        junkname = DebugDir + "/" + junkname;

        if (ShowMF_Verbose>0 && ParallelDescriptor::IOProcessor()) {
            cout << "   ******************************  Debug: writing " << junkname << endl;
        }

        FArrayBox sub(box,mf.nComp());

        MultiFab mfg(BoxArray(mf.boxArray()).grow(mf.nGrow()),mf.nComp(),0);
        for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        {
            mfg[mfi].copy(mf[mfi]);
        }
        mfg.copy(sub,0,0,mf.nComp());

        if (ShowMF_Check_Nans)
        {
            BL_ASSERT(!sub.contains_nan(box,0,mf.nComp()));
        }
        std::ofstream os;
        os.open(junkname.c_str());
        sub.writeOn(os);
        os.close();
        FArrayBox::setFormat(saved_format);
    }
}

void
showMF(const std::string&   mySet,
       const MultiFab&      mf,
       const std::string&   name,
       int                  lev = -1,
       int                  iter = -1) // Default value = no append 2nd integer
{
    if (ShowMF_Sets.count(mySet)>0)
    {
        const FABio::Format saved_format = FArrayBox::getFormat();
        FArrayBox::setFormat(ShowMF_Fab_Format);
        std::string DebugDir(ShowMF_Dir);
        if (ParallelDescriptor::IOProcessor())
            if (!BoxLib::UtilCreateDirectory(DebugDir, 0755))
                BoxLib::CreateDirectoryFailed(DebugDir);
        ParallelDescriptor::Barrier();

        std::string junkname = name;
        if (lev>=0) {
            junkname = BoxLib::Concatenate(junkname+"_",lev,1);
        }
        if (iter>=0) {
            junkname = BoxLib::Concatenate(junkname+"_",iter,1);
        }
        junkname = DebugDir + "/" + junkname;

        if (ShowMF_Verbose>0 && ParallelDescriptor::IOProcessor()) {
            cout << "   ******************************  Debug: writing " << junkname << endl;
        }

        if (ShowMF_Check_Nans)
        {
            for (MFIter mfi(mf); mfi.isValid(); ++mfi)
            {
                BL_ASSERT(!mf[mfi].contains_nan(mfi.validbox(),0,mf.nComp()));
            }
        }
        VisMF::Write(mf,junkname);
        FArrayBox::setFormat(saved_format);
    }
}

void
showMCDD(const std::string&   mySet,
         const DDOp&          mcddop,
         const std::string&   name)
{
    if (ShowMF_Sets.count(mySet)>0)
    {
        std::string DebugDir(ShowMF_Dir);
        if (ParallelDescriptor::IOProcessor())
            if (!BoxLib::UtilCreateDirectory(DebugDir, 0755))
                BoxLib::CreateDirectoryFailed(DebugDir);
        ParallelDescriptor::Barrier();

        std::string junkname = DebugDir + "/" + name;

        if (ShowMF_Verbose>0 && ParallelDescriptor::IOProcessor()) {
            cout << "   ******************************  Debug: writing " << junkname << endl;
        }
        mcddop.Write(junkname);
    }
}

HeatTransfer::FPLoc 
HeatTransfer::fpi_phys_loc (int p_bc)
{
    //
    // Location of data that FillPatchIterator returns at physical boundaries
    //
    if (p_bc == EXT_DIR || p_bc == HOEXTRAP || p_bc == FOEXTRAP)
    {
        return HT_Edge;
    }
    return HT_Center;
}
    
void
HeatTransfer::center_to_edge_fancy (const FArrayBox& cfab,
                                    FArrayBox&       efab,
                                    const Box&       ccBox,
                                    int              sComp,
                                    int              dComp,
                                    int              nComp,
                                    const Box&       domain,
                                    const FPLoc&     bc_lo,
                                    const FPLoc&     bc_hi)
{
    const Box&      ebox = efab.box();
    const IndexType ixt  = ebox.ixType();

    BL_ASSERT(!(ixt.cellCentered()) && !(ixt.nodeCentered()));

    int dir = -1;
    for (int d = 0; d < BL_SPACEDIM; d++)
        if (ixt.test(d))
            dir = d;

    BL_ASSERT(BoxLib::grow(ccBox,-BoxLib::BASISV(dir)).contains(BoxLib::enclosedCells(ebox)));
    BL_ASSERT(sComp+nComp <= cfab.nComp() && dComp+nComp <= efab.nComp());
    //
    // Exclude unnecessary cc->ec calcs
    //
    Box ccVBox = ccBox;
    if (bc_lo!=HT_Center)
        ccVBox.setSmall(dir,std::max(domain.smallEnd(dir),ccVBox.smallEnd(dir)));
    if (bc_hi!=HT_Center)
        ccVBox.setBig(dir,std::min(domain.bigEnd(dir),ccVBox.bigEnd(dir)));
    //
    // Shift cell-centered data to edges
    //
    const int isharm = def_harm_avg_cen2edge?1:0;
    FORT_CEN2EDG(ccVBox.loVect(),ccVBox.hiVect(),
                 ARLIM(cfab.loVect()),ARLIM(cfab.hiVect()),cfab.dataPtr(sComp),
                 ARLIM(efab.loVect()),ARLIM(efab.hiVect()),efab.dataPtr(dComp),
                 &nComp, &dir, &isharm);
    //
    // Fix boundary...i.e. fill-patched data in cfab REALLY lives on edges
    //
    if ( !(domain.contains(ccBox)) )
    {
        if (bc_lo==HT_Edge)
        {
            BoxList gCells = BoxLib::boxDiff(ccBox,domain);
            if (gCells.ok())
            {
                const int inc = +1;
                FArrayBox ovlpFab;
                for (BoxList::const_iterator bli = gCells.begin(), end = gCells.end();
                     bli != end;
                     ++bli)
                {
                    if (bc_lo == HT_Edge)
                    {
                        ovlpFab.resize(*bli,nComp);
                        ovlpFab.copy(cfab,sComp,0,nComp);
                        ovlpFab.shiftHalf(dir,inc);
                        efab.copy(ovlpFab,0,dComp,nComp);
                    }
                }
            }
        }
        if (bc_hi==HT_Edge)
        {
            BoxList gCells = BoxLib::boxDiff(ccBox,domain);
            if (gCells.ok())
            {
                const int inc = -1;
                FArrayBox ovlpFab;
                for (BoxList::const_iterator bli = gCells.begin(), end = gCells.end();
                     bli != end;
                     ++bli)
                {
                    if (bc_hi == HT_Edge)
                    {
                        ovlpFab.resize(*bli,nComp);
                        ovlpFab.copy(cfab,sComp,0,nComp);
                        ovlpFab.shiftHalf(dir,inc);
                        efab.copy(ovlpFab,0,dComp,nComp);
                    }
                }
            }
        }
    }
}    

void
HeatTransfer::variableCleanUp ()
{
    NavierStokes::variableCleanUp();

    delete chemSolve;
    chemSolve = 0;

    mcdd_nu1.clear();
    mcdd_nu2.clear();
    ShowMF_Sets.clear();
    auxDiag_names.clear();
    typical_values.clear();

#ifdef PARTICLES
    delete HTPC;
    HTPC = 0;
#endif
}

HeatTransfer::HeatTransfer ()
{
    if (!init_once_done)
        init_once();

    if (!do_temp)
        BoxLib::Abort("do_temp MUST be true");

    if (!have_divu)
        BoxLib::Abort("have_divu MUST be true");

    if (!have_dsdt)
        BoxLib::Abort("have_dsdt MUST be true");

    EdgeState              = 0;
    EdgeFlux               = 0;
    SpecDiffusionFluxn     = 0;
    SpecDiffusionFluxnp1   = 0;
    FillPatchedOldState_ok = true;
}

HeatTransfer::HeatTransfer (Amr&            papa,
                            int             lev,
                            const Geometry& level_geom,
                            const BoxArray& bl,
                            Real            time)
    :
    NavierStokes(papa,lev,level_geom,bl,time),
    //
    // Make room for all components except velocities in aux_boundary_data_old.
    //
    aux_boundary_data_old(bl,HYP_GROW,desc_lst[State_Type].nComp()-BL_SPACEDIM,level_geom),
    //
    // Only save Density & RhoH in aux_boundary_data_new in components 0 & 1.
    //
    aux_boundary_data_new(bl,LinOp_grow,2,level_geom),
    FillPatchedOldState_ok(true)
{
    if (!init_once_done)
        init_once();

    if (!do_temp)
        BoxLib::Abort("do_temp MUST be true");

    if (!have_divu)
        BoxLib::Abort("have_divu MUST be true");

    if (!have_dsdt)
        BoxLib::Abort("have_dsdt MUST be true");

    const int nGrow       = 0;
    const int nEdgeStates = desc_lst[State_Type].nComp();
    diffusion->allocFluxBoxesLevel(EdgeState,nGrow,nEdgeStates);
    diffusion->allocFluxBoxesLevel(EdgeFlux,nGrow,nEdgeStates);
    if (nspecies>0 && !unity_Le)
    {
	diffusion->allocFluxBoxesLevel(SpecDiffusionFluxn,  nGrow,nspecies+3);
	diffusion->allocFluxBoxesLevel(SpecDiffusionFluxnp1,nGrow,nspecies+3);
        sumSpecFluxDotGradHn.define(grids,1,0,Fab_allocate);
        sumSpecFluxDotGradHnp1.define(grids,1,0,Fab_allocate);
    }

    for (std::map<std::string,Array<std::string> >::iterator it = auxDiag_names.begin(), end = auxDiag_names.end();
         it != end; ++it)
    {
        auxDiag[it->first] = new MultiFab(grids,it->second.size(),0);
        auxDiag[it->first]->setVal(0);
    }

    if (do_mcdd)
        MCDDOp.define(grids,Domain(),crse_ratio);

    Forcing.define(grids,nspecies+1,nGrowAdvForcing,Fab_allocate);
    Dn.define(grids,nspecies+2,nGrowAdvForcing,Fab_allocate);
    Dhat.define(grids,nspecies+2,nGrowAdvForcing,Fab_allocate);
    Dnp1.define(grids,nspecies+2,nGrowAdvForcing,Fab_allocate);
    DDn.define(grids,1,nGrowAdvForcing,Fab_allocate);
    DDnp1.define(grids,1,nGrowAdvForcing,Fab_allocate);

    // HACK for debugging
    if (level==0)
        stripBox = getStrip(geom);

}

HeatTransfer::~HeatTransfer ()
{
    diffusion->removeFluxBoxesLevel(EdgeState);
    diffusion->removeFluxBoxesLevel(EdgeFlux);
    if (nspecies>0 && !unity_Le)    
    {
	diffusion->removeFluxBoxesLevel(SpecDiffusionFluxn);    
	diffusion->removeFluxBoxesLevel(SpecDiffusionFluxnp1);    
        sumSpecFluxDotGradHn.clear();
        sumSpecFluxDotGradHnp1.clear();
    }

    for (std::map<std::string,MultiFab*>::iterator it = auxDiag.begin(), end = auxDiag.end();
         it != end;
         ++it)
    {
        delete it->second;
    }
}

void
HeatTransfer::init_once ()
{
    //
    // Computes the static variables unique to HeatTransfer.
    // Check that (some) things are set up correctly.
    //
    int dummy_State_Type;

    int have_temp = isStateVariable("temp", dummy_State_Type, Temp);

    have_temp = have_temp && State_Type == dummy_State_Type;
    have_temp = have_temp && isStateVariable("rhoh", dummy_State_Type, RhoH);
    have_temp = have_temp && State_Type == dummy_State_Type;

    have_trac = isStateVariable("tracer", dummy_State_Type, Trac);
    have_trac = have_trac && State_Type == dummy_State_Type;

    have_rhort = isStateVariable("RhoRT", dummy_State_Type, RhoRT);
    have_rhort = have_rhort && State_Type == dummy_State_Type;

    if (!have_temp)
        BoxLib::Abort("HeatTransfer::init_once(): RhoH & Temp must both be the state");
    
    if (!have_rhort && verbose && ParallelDescriptor::IOProcessor())
        BoxLib::Warning("HeatTransfer::init_once(): RhoRT being stored in the Tracer slot");
    
    if (Temp < RhoH)
        BoxLib::Abort("HeatTransfer::init_once(): must have RhoH < Temp");
    //
    // Temperature must be non-conservative, rho*h must be conservative.
    //
    if (advectionType[Temp] == Conservative)
        BoxLib::Abort("HeatTransfer::init_once(): Temp must be non-conservative");

    if (advectionType[RhoH] != Conservative)
        BoxLib::Abort("HeatTransfer::init_once(): RhoH must be conservative");
    //
    // Species checks.
    //
    BL_ASSERT(Temp > RhoH && RhoH > Density);
    //
    // Here we want to count relative to Density instead of relative
    // to RhoH, so we can put RhoH after the species instead of before.  
    // This logic should work in both cases.
    //
    first_spec =  Density + 1;
    last_spec  = first_spec + getChemSolve().numSpecies() - 1;
    
    for (int i = first_spec; i <= last_spec; i++)
        if (advectionType[i] != Conservative)
            BoxLib::Error("HeatTransfer::init_once: species must be conservative");
    
    int diffuse_spec = is_diffusive[first_spec];
    for (int i = first_spec+1; i <= last_spec; i++)
        if (is_diffusive[i] != diffuse_spec)
            BoxLib::Error("HeatTransfer::init_once: Le != 1; diffuse");
    //
    // Load integer pointers into Fortran common, reqd for proper ICs.
    //
    const int density = (int)Density;

    FORT_SET_SCAL_NUMB(&density, &Temp, &Trac, &RhoH, &first_spec, &last_spec);
    //
    // Load constants into Fortran common to compute viscosities, etc.
    //
    const int var_visc = (constant_mu_val     == -1 ? 1 : 0);
    const int var_cond = (constant_lambda_val == -1 ? 1 : 0);
    const int var_diff = (constant_rhoD_val   == -1 ? 1 : 0);
    
    FORT_SET_HT_VISC_COMMON(&var_visc, &constant_mu_val,
                            &var_cond, &constant_lambda_val,
                            &var_diff, &constant_rhoD_val,
                            &prandtl,  &schmidt, &unity_Le);

    //
    // make space for typical values
    //
    typical_values.resize(NUM_STATE,-1); // -ve means don't use for anything
    //
    // Get universal gas constant from Fortran.
    //
    rgas = getChemSolve().getRuniversal();
    P1atm_MKS = getChemSolve().getP1atm_MKS();

    if (rgas <= 0.0)
    {
        std::cerr << "HeatTransfer::init_once(): bad rgas: " << rgas << '\n';
        BoxLib::Abort();
    }
    if (P1atm_MKS <= 0.0)
    {
        std::cerr << "HeatTransfer::init_once(): bad P1atm_MKS: " << P1atm_MKS << '\n';
        BoxLib::Abort();
    }
    //
    // Chemistry.
    //
    int ydot_good = RhoYdot_Type >= 0 && RhoYdot_Type <desc_lst.size()
        && RhoYdot_Type != Divu_Type
        && RhoYdot_Type != Dsdt_Type
        && RhoYdot_Type != State_Type;
    
    if (!ydot_good)
        BoxLib::Error("HeatTransfer::init_once(): need RhoYdot_Type if do_chemistry");
    
    const StateDescriptor& ydot_cell = desc_lst[RhoYdot_Type];
    int nydot = ydot_cell.nComp();
    if (nydot < nspecies)
        BoxLib::Error("HeatTransfer::init_once(): RhoYdot_Type needs nspecies components");
    //
    // Enforce Le = 1, unless !unity_Le
    //
    if (unity_Le && (schmidt != prandtl) )
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "**************** WARNING ***************\n"
                      << "HeatTransfer::init_once() : currently must have"
                      << "equal Schmidt and Prandtl numbers unless !unity_Le.\n"
                      << "Setting Schmidt = Prandtl\n"
                      << "**************** WARNING ***************\n";
    
        schmidt = prandtl;
    }
    //
    // We are done.
    //
    num_state_type = desc_lst.size();

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::init_once(): num_state_type = " << num_state_type << '\n';

    ParmParse pp("ht");

    pp.query("plot_reactions",plot_reactions);
    if (plot_reactions)
    {
        auxDiag_names["REACTIONS"].resize(getChemSolve().numReactions());
        for (int i = 0; i < auxDiag_names["REACTIONS"].size(); ++i)
            auxDiag_names["REACTIONS"][i] = BoxLib::Concatenate("R",i+1);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "***** Make sure to increase amr.regrid_int !!!!!" << '\n';
    }

    pp.query("plot_consumption",plot_consumption);
    pp.query("plot_auxDiags",plot_consumption); // This is for backward comptibility - FIXME
    if (plot_consumption)
    {
        auxDiag_names["CONSUMPTION"].resize(consumptionName.size());
        for (int j=0; j<consumptionName.size(); ++j)
        {
            auxDiag_names["CONSUMPTION"][j] = consumptionName[j] + "_ConsumptionRate";
        }
    }

    pp.query("plot_heat_release",plot_heat_release);
    if (plot_heat_release)
    {
        auxDiag_names["HEATRELEASE"].resize(1);
        auxDiag_names["HEATRELEASE"][0] = "HeatRelease";
    }
    
    pp.query("new_T_threshold",new_T_threshold);

    init_once_done = 1;
}

void
HeatTransfer::restart (Amr&          papa,
                       std::istream& is,
                       bool          bReadSpecial)
{

    NavierStokes::restart(papa,is,bReadSpecial);
    //
    // Make room for all components except velocities in aux_boundary_data_old.
    //
    aux_boundary_data_old.initialize(grids,HYP_GROW,desc_lst[State_Type].nComp()-BL_SPACEDIM,Geom());
    //
    // Only save Density & RhoH in aux_boundary_data_new in components 0 & 1.
    //
    aux_boundary_data_new.initialize(grids,LinOp_grow,2,Geom());

    FillPatchedOldState_ok = true;
    
    // set_overdetermined_boundary_cells(state[State_Type].curTime());

    BL_ASSERT(EdgeState == 0);
    BL_ASSERT(EdgeFlux == 0);
    const int nGrow       = 0;
    const int nEdgeStates = desc_lst[State_Type].nComp();
    diffusion->allocFluxBoxesLevel(EdgeState,nGrow,nEdgeStates);
    diffusion->allocFluxBoxesLevel(EdgeFlux,nGrow,nEdgeStates);
    
    if (nspecies>0 && !unity_Le)
    {
	BL_ASSERT(SpecDiffusionFluxn == 0);
	BL_ASSERT(SpecDiffusionFluxnp1 == 0);
	diffusion->allocFluxBoxesLevel(SpecDiffusionFluxn,  nGrow,nspecies+3);
	diffusion->allocFluxBoxesLevel(SpecDiffusionFluxnp1,nGrow,nspecies+3);
        sumSpecFluxDotGradHn.define(grids,1,0,Fab_allocate);
        sumSpecFluxDotGradHnp1.define(grids,1,0,Fab_allocate);
    }

    for (std::map<std::string,Array<std::string> >::iterator it = auxDiag_names.begin(), end = auxDiag_names.end();
         it != end; ++it)
    {
        auxDiag[it->first] = new MultiFab(grids,it->second.size(),0);
        auxDiag[it->first]->setVal(0);
    }

    if (do_mcdd)
        MCDDOp.define(grids,Domain(),crse_ratio);

    Forcing.define(grids,nspecies+1,nGrowAdvForcing,Fab_allocate);
    Dn.define(grids,nspecies+2,nGrowAdvForcing,Fab_allocate);
    Dhat.define(grids,nspecies+2,nGrowAdvForcing,Fab_allocate);
    Dnp1.define(grids,nspecies+2,nGrowAdvForcing,Fab_allocate);
    DDn.define(grids,1,nGrowAdvForcing,Fab_allocate);
    DDnp1.define(grids,1,nGrowAdvForcing,Fab_allocate);

    // HACK for debugging
    if (level==0)
        stripBox = getStrip(geom);

    // Deal with typical values
    set_typical_values(true);
}

void
HeatTransfer::set_typical_values(bool restart)
{
    if (level==0) {

        int nComp = typical_values.size();

        if (restart) {

            BL_ASSERT(nComp==NUM_STATE);
            for (int i=0; i<nComp; ++i) {
                typical_values[i] = -1;
            }
            
            if (ParallelDescriptor::IOProcessor())
            {
                
                const std::string tvfile = parent->theRestartFile() + "/" + typical_values_filename;
                std::ifstream tvis;
                tvis.open(tvfile.c_str(),std::ios::in);
                
                if (tvis.good()) {
                    FArrayBox tvfab;
                    tvfab.readFrom(tvis);
                    if (tvfab.nComp() != typical_values.size()) {
                        BoxLib::Abort("Typical values file has wrong number of components");
                    }
                    for (int i=0; i<typical_values.size(); ++i) {
                        typical_values[i] = tvfab.dataPtr()[i];
                    }
                }
            }

            ParallelDescriptor::ReduceRealMax(typical_values.dataPtr(),nComp); //FIXME: better way?
        }

        // Check fortan common values, override values set above if fortran values > 0
        Array<Real> tvTmp(nComp,0);
        FORT_GETTYPICALVALS(tvTmp.dataPtr(), &nComp);
        ParallelDescriptor::ReduceRealMax(tvTmp.dataPtr(),nComp);
        
        for (int i=0; i<nComp; ++i) {
            if (tvTmp[i]>0) {
                typical_values[i] = tvTmp[i];
            }
        }
        FORT_SETTYPICALVALS(typical_values.dataPtr(), &nComp);
	ParallelDescriptor::ReduceRealMax(typical_values.dataPtr(),nComp);

        if (ParallelDescriptor::IOProcessor())
        {
            cout << "Typical vals: " << endl;
            cout << "\tVelocity: ";
            for (int i=0; i<BL_SPACEDIM; ++i) {
                cout << typical_values[i] << " ";
            }
            cout << endl;
            cout << "\tDensity: " << typical_values[Density] << endl;
            cout << "\tTemp: "    << typical_values[Temp]    << endl;
            cout << "\tRhoH: "    << typical_values[RhoH]    << endl;
            const Array<std::string>& names = getChemSolve().speciesNames();
            for (int i=0; i<nspecies; ++i) {
                cout << "\tY_" << names[i] << ": " << typical_values[first_spec+i] << endl;
            }
        }
            
        // Verify good values for Density, Temp, RhoH, Y -- currently only needed for mcdd problems
        if (do_mcdd && ParallelDescriptor::IOProcessor()) {
            for (int i=BL_SPACEDIM; i<nComp; ++i) {
                if (i!=Trac && i!=RhoRT && typical_values[i]<=0) {
                    cout << "component: " << i << " of " << nComp << endl;
                    BoxLib::Abort("Must have non-zero typical values");
                }
            }
        }
    }
}

Real
HeatTransfer::estTimeStep ()
{
    Real estdt = NavierStokes::estTimeStep();

    if (fixed_dt > 0.0 || !divu_ceiling)
        //
        // The n-s function did the right thing in this case.
        //
        return estdt;

    Real dt, ns_estdt = estdt, divu_dt = 1.0e20;

    const int   n_grow   = 1;
    const Real  cur_time = state[State_Type].curTime();
    const Real* dx       = geom.CellSize();
    MultiFab*   dsdt     = getDsdt(0,cur_time);
    MultiFab*   divu     = getDivCond(0,cur_time);

    FArrayBox area[BL_SPACEDIM], volume;

    for (FillPatchIterator U_fpi(*this,*divu,n_grow,cur_time,State_Type,Xvel,BL_SPACEDIM);
         U_fpi.isValid();
         ++U_fpi)
    {
        const int        i   = U_fpi.index();
        FArrayBox&       U   = U_fpi();
        const FArrayBox& Rho = (*rho_ctime)[i];
        const int*       lo  = grids[i].loVect();
        const int*       hi  = grids[i].hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            geom.GetFaceArea(area[dir],grids,i,dir,GEOM_GROW);

        geom.GetVolume(volume, grids, i, GEOM_GROW);

        DEF_CLIMITS((*divu)[i],sdat,slo,shi);
        DEF_CLIMITS(Rho,rhodat,rholo,rhohi);
        DEF_CLIMITS(U,vel,ulo,uhi);

        DEF_CLIMITS(volume,vol,v_lo,v_hi);

        DEF_CLIMITS(area[0],areax,ax_lo,ax_hi);
        DEF_CLIMITS(area[1],areay,ay_lo,ay_hi);
#if (BL_SPACEDIM==3) 
        DEF_CLIMITS(area[2],areaz,az_lo,az_hi)
#endif
        FORT_EST_DIVU_DT(divu_ceiling,&divu_dt_factor,
                         dx,sdat,ARLIM(slo),ARLIM(shi),
                         (*dsdt)[i].dataPtr(),
                         rhodat,ARLIM(rholo),ARLIM(rhohi),
                         vel,ARLIM(ulo),ARLIM(uhi),
                         vol,ARLIM(v_lo),ARLIM(v_hi),
                         areax,ARLIM(ax_lo),ARLIM(ax_hi),
                         areay,ARLIM(ay_lo),ARLIM(ay_hi),
#if (BL_SPACEDIM==3) 
                         areaz,ARLIM(az_lo),ARLIM(az_hi),
#endif 
                         lo,hi,&dt,&min_rho_divu_ceiling);

        divu_dt = std::min(divu_dt,dt);
    }

    delete divu;
    delete dsdt;

    ParallelDescriptor::ReduceRealMin(divu_dt);

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "HeatTransfer::estTimeStep(): estdt, divu_dt = " 
                  << estdt << ", " << divu_dt << '\n';
    }

    estdt = std::min(estdt, divu_dt);

    if (estdt < ns_estdt && verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::estTimeStep(): timestep reduced from " 
                  << ns_estdt << " to " << estdt << '\n';

    return estdt;
}

void
HeatTransfer::checkTimeStep (Real dt)
{
    if (fixed_dt > 0.0 || !divu_ceiling) 
        return;

    const int   n_grow    = 1;
    const Real  cur_time  = state[State_Type].curTime();
    const Real* dx        = geom.CellSize();
    MultiFab*   dsdt      = getDsdt(0,cur_time);
    MultiFab*   divu      = getDivCond(0,cur_time);

    FArrayBox area[BL_SPACEDIM], volume;

    for (FillPatchIterator U_fpi(*this,*divu,n_grow,cur_time,State_Type,Xvel,BL_SPACEDIM);
         U_fpi.isValid();
         ++U_fpi)
    {
        const int        i   = U_fpi.index();
        FArrayBox&       U   = U_fpi();
        const FArrayBox& Rho = (*rho_ctime)[i];
        const int*       lo  = grids[i].loVect();
        const int*       hi  = grids[i].hiVect();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            geom.GetFaceArea(area[dir],grids,i,dir,GEOM_GROW);

        geom.GetVolume(volume, grids, i, GEOM_GROW);

        DEF_LIMITS((*divu)[U_fpi],sdat,slo,shi);
        DEF_CLIMITS(Rho,rhodat,rholo,rhohi);
        DEF_CLIMITS(U,vel,ulo,uhi);

        DEF_CLIMITS(volume,vol,v_lo,v_hi);
        DEF_CLIMITS(area[0],areax,ax_lo,ax_hi);
        DEF_CLIMITS(area[1],areay,ay_lo,ay_hi);

#if (BL_SPACEDIM==3) 
        DEF_CLIMITS(area[2],areaz,az_lo,az_hi);
#endif
        FORT_CHECK_DIVU_DT(divu_ceiling,&divu_dt_factor,
                           dx,sdat,ARLIM(slo),ARLIM(shi),
			   (*dsdt)[U_fpi].dataPtr(),
                           rhodat,ARLIM(rholo),ARLIM(rhohi),
                           vel,ARLIM(ulo),ARLIM(uhi),
                           vol,ARLIM(v_lo),ARLIM(v_hi),
                           areax,ARLIM(ax_lo),ARLIM(ax_hi),
                           areay,ARLIM(ay_lo),ARLIM(ay_hi),
#if (BL_SPACEDIM==3) 
                           areaz,ARLIM(az_lo),ARLIM(az_hi),
#endif 
                           lo,hi,&dt,&min_rho_divu_ceiling);
    }

    delete dsdt;
    delete divu;
}

void
HeatTransfer::setTimeLevel (Real time,
                            Real dt_old,
                            Real dt_new)
{
    NavierStokes::setTimeLevel(time, dt_old, dt_new);    

    state[RhoYdot_Type].setTimeLevel(time,dt_old,dt_new);

    state[FuncCount_Type].setTimeLevel(time,dt_old,dt_new);
}

//
// This (minus the NEWMECH stuff) is copied from NavierStokes.cpp
//

void
HeatTransfer::initData ()
{
    //
    // Initialize the state and the pressure.
    //
    int         ns       = NUM_STATE - BL_SPACEDIM;
    const Real* dx       = geom.CellSize();
    MultiFab&   S_new    = get_new_data(State_Type);
    MultiFab&   P_new    = get_new_data(Press_Type);
    const Real  cur_time = state[State_Type].curTime();

#ifdef BL_USE_NEWMECH
    //
    // This code has a few drawbacks.  It assumes that the physical
    // domain size of the current problem is the same as that of the
    // one that generated the pltfile.  It also assumes that the pltfile
    // has at least as many levels as does the current problem.  If
    // either of these are false this code is likely to core dump.
    //
    ParmParse pp("ht");

    std::string pltfile;
    pp.query("pltfile", pltfile);
    if (pltfile.empty())
        BoxLib::Abort("You must specify `pltfile'");
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "initData: reading data from: " << pltfile << '\n';

    DataServices::SetBatchMode();
    FileType fileType(NEWPLT);
    DataServices dataServices(pltfile, fileType);

    if (!dataServices.AmrDataOk())
        //
        // This calls ParallelDescriptor::EndParallel() and exit()
        //
        DataServices::Dispatch(DataServices::ExitRequest, NULL);
    
    AmrData&                  amrData     = dataServices.AmrDataRef();
    const int                 nspecies    = getChemSolve().numSpecies();
    const Array<std::string>& names       = getChemSolve().speciesNames();   
    Array<std::string>        plotnames   = amrData.PlotVarNames();

    int idT = -1, idX = -1;
    for (int i = 0; i < plotnames.size(); ++i)
    {
        if (plotnames[i] == "temp")       idT = i;
        if (plotnames[i] == "x_velocity") idX = i;
    }
    //
    // In the plotfile the mass fractions directly follow the velocities.
    //
    int idSpec = idX + BL_SPACEDIM;

    for (int i = 0; i < BL_SPACEDIM; i++)
    {
        amrData.FillVar(S_new, level, plotnames[idX+i], Xvel+i);
        amrData.FlushGrids(idX+i);
    }
    amrData.FillVar(S_new, level, plotnames[idT], Temp);
    amrData.FlushGrids(idT);

    for (int i = 0; i < nspecies; i++)
    {
        amrData.FillVar(S_new, level, plotnames[idSpec+i], first_spec+i);
        amrData.FlushGrids(idSpec+i);
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "initData: finished init from pltfile" << '\n';
#endif

    for (MFIter snewmfi(S_new); snewmfi.isValid(); ++snewmfi)
    {
        BL_ASSERT(grids[snewmfi.index()] == snewmfi.validbox());

        P_new[snewmfi].setVal(0);

        const int  i       = snewmfi.index();
        RealBox    gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
        const int* lo      = snewmfi.validbox().loVect();
        const int* hi      = snewmfi.validbox().hiVect();
        const int* s_lo    = S_new[snewmfi].loVect();
        const int* s_hi    = S_new[snewmfi].hiVect();
        const int* p_lo    = P_new[snewmfi].loVect();
        const int* p_hi    = P_new[snewmfi].hiVect();

#ifdef BL_USE_NEWMECH
        FORT_INITDATANEWMECH (&level,&cur_time,lo,hi,&ns,
                              S_new[snewmfi].dataPtr(Xvel),
                              S_new[snewmfi].dataPtr(BL_SPACEDIM),
                              ARLIM(s_lo), ARLIM(s_hi),
                              P_new[snewmfi].dataPtr(),
                              ARLIM(p_lo), ARLIM(p_hi),
                              dx,gridloc.lo(),gridloc.hi() );
#else
        FORT_INITDATA (&level,&cur_time,lo,hi,&ns,
                       S_new[snewmfi].dataPtr(Xvel),
                       S_new[snewmfi].dataPtr(BL_SPACEDIM),
                       ARLIM(s_lo), ARLIM(s_hi),
                       P_new[snewmfi].dataPtr(),
                       ARLIM(p_lo), ARLIM(p_hi),
                       dx,gridloc.lo(),gridloc.hi() );
#endif
    }

    showMFsub("1D",S_new,stripBox,"1D_S",level);

#ifdef BL_USE_VELOCITY
    //
    // We want to add the velocity from the supplied plotfile
    // to what we already put into the velocity field via FORT_INITDATA.
    //
    // This code has a few drawbacks.  It assumes that the physical
    // domain size of the current problem is the same as that of the
    // one that generated the pltfile.  It also assumes that the pltfile
    // has at least as many levels (with the same refinement ratios) as does
    // the current problem.  If either of these are false this code is
    // likely to core dump.
    //
    ParmParse pp("ht");

    std::string velocity_plotfile;
    pp.query("velocity_plotfile", velocity_plotfile);

    if (!velocity_plotfile.empty())
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "initData: reading data from: " << velocity_plotfile << '\n';

        DataServices::SetBatchMode();
        Amrvis::FileType fileType(Amrvis::NEWPLT);
        DataServices dataServices(velocity_plotfile, fileType);

        if (!dataServices.AmrDataOk())
            //
            // This calls ParallelDescriptor::EndParallel() and exit()
            //
            DataServices::Dispatch(DataServices::ExitRequest, NULL);

        AmrData&                  amrData   = dataServices.AmrDataRef();
        //const int                 nspecies  = getChemSolve().numSpecies();
        //const Array<std::string>& names     = getChemSolve().speciesNames();   
        Array<std::string>        plotnames = amrData.PlotVarNames();

        if (amrData.FinestLevel() < level)
            BoxLib::Abort("initData: not enough levels in plotfile");

        if (amrData.ProbDomain()[level] != Domain())
            BoxLib::Abort("initData: problem domains do not match");
    
        int idX = -1;
        for (int i = 0; i < plotnames.size(); ++i)
            if (plotnames[i] == "x_velocity") idX = i;

        if (idX == -1)
            BoxLib::Abort("Could not find velocity fields in supplied velocity_plotfile");

        MultiFab tmp(S_new.boxArray(), 1, 0);
        for (int i = 0; i < BL_SPACEDIM; i++)
        {
            amrData.FillVar(tmp, level, plotnames[idX+i], 0);
            for (MFIter mfi(tmp); mfi.isValid(); ++mfi)
                S_new[mfi].plus(tmp[mfi], tmp[mfi].box(), 0, Xvel+i, 1);
            amrData.FlushGrids(idX+i);
        }

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "initData: finished init from velocity_plotfile" << '\n';
    }
#endif /*BL_USE_VELOCITY*/

    make_rho_prev_time();
    make_rho_curr_time();
    //
    // Initialize other types.
    //
    initDataOtherTypes();
    //
    // Initialize divU and dSdt.
    //
    if (have_divu)
    {
        const Real dt       = 1.0;
        const Real dtin     = -1.0; // Dummy value denotes initialization.
        const Real cur_time = state[Divu_Type].curTime();
        MultiFab&  Divu_new = get_new_data(Divu_Type);

        state[State_Type].setTimeLevel(cur_time,dt,dt);

        calcDiffusivity(cur_time);
        calc_divu(cur_time,dtin,Divu_new);

        if (have_dsdt)
            get_new_data(Dsdt_Type).setVal(0);
    }

    if (state[Press_Type].descriptor()->timeType() == StateDescriptor::Point) 
    {
        get_new_data(Dpdt_Type).setVal(0);
    }

    is_first_step_after_regrid = false;
    old_intersect_new          = grids;

#ifdef PARTICLES
    if (level == 0)
    {
        if (HTPC == 0)
        {
            HTPC = new HTParticleContainer(parent);
            //
            // Make sure to call RemoveParticles() on exit.
            //
            BoxLib::ExecOnFinalize(RemoveParticles);
        }

        HTPC->SetVerbose(pverbose);

        if (!particle_init_file.empty())
        {
            HTPC->InitFromAsciiFile(particle_init_file,0);
        }
    }
#endif
}

void
HeatTransfer::initDataOtherTypes ()
{
    // to be consistent with Strang code, compute enthalpy with the EOS
    const Real cur_time  = state[State_Type].curTime();
    {
        MultiFab rhoh(grids,1,0);
	compute_rhohmix(cur_time,rhoh);
        get_new_data(State_Type).copy(rhoh,0,RhoH,1);
    }
    //
    // Assume that by now, S_new has "good" data
    //
    MultiFab& R = get_new_data(RhoYdot_Type);

    // It is our current belief that at the beginning of a simulation we want omegadot=0
    // This is how the Strang code does it
    //    compute_instantaneous_reaction_rates(R,get_new_data(State_Type),nGrow);
    get_new_data(RhoYdot_Type).setVal(0);

    showMFsub("1D",R,stripBox,"1D_dd_idot_R",level);

    get_new_data(FuncCount_Type).setVal(1);

    setThermoPress(state[State_Type].curTime());
}

void
HeatTransfer::compute_instantaneous_reaction_rates(MultiFab&       R,
                                                   const MultiFab& S,
                                                   int             nGrow,
                                                   HowToFillGrow   how)
{
    Real p_amb, dpdt_factor;
    FORT_GETPAMB(&p_amb, &dpdt_factor);
    const Real Patm = p_amb / P1atm_MKS;

    BL_ASSERT((nGrow==0)  ||  (how == HT_ZERO_GROW_CELLS) || (how == HT_EXTRAP_GROW_CELLS));

    if ((nGrow>0) && (how == HT_ZERO_GROW_CELLS))
    {
        R.setBndry(0,0,nspecies);
    }
    
    int sCompRhoH = RhoH;
    int sCompRhoY = first_spec;
    int sCompT    = Temp;
    int sCompRhoYdot = 0;
    
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const FArrayBox& rhoY = S[mfi];
        const FArrayBox& rhoH = S[mfi];
        const FArrayBox& T    = S[mfi];
        const Box& box = mfi.validbox();
        FArrayBox& rhoYdot = R[mfi];

        getChemSolve().reactionRateRhoY(rhoYdot,rhoY,rhoH,T,Patm,box,sCompRhoY,sCompRhoH,sCompT,sCompRhoYdot);
    }
    
    
    if ((nGrow>0) && (how == HT_EXTRAP_GROW_CELLS))
    {
        const int N = R.IndexMap().size();
        
#ifdef BL_USE_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N; i++)
        {
            const int  k   = R.IndexMap()[i];
            FArrayBox& r   = R[k];
            const Box& box = R.box(k);
            FORT_VISCEXTRAP(r.dataPtr(),ARLIM(r.loVect()),ARLIM(r.hiVect()),
                            box.loVect(),box.hiVect(),&nspecies);
        }
        R.FillBoundary(0,nspecies);
        //
        // Note: this is a special periodic fill in that we want to
        // preserve the extrapolated grow values when periodic --
        // usually we preserve only valid data.  The scheme relies on
        // the fact that there is good data in the "non-periodic" grow cells.
        // ("good" data produced via VISCEXTRAP above)
        //
        geom.FillPeriodicBoundary(R,0,nspecies,true);
    }
} 

void
HeatTransfer::init (AmrLevel& old)
{
    NavierStokes::init(old);

    HeatTransfer* oldht    = (HeatTransfer*) &old;
    const Real    cur_time = oldht->state[State_Type].curTime();
    //
    // Get best ydot data.
    //
    MultiFab& Ydot = get_new_data(RhoYdot_Type);

    for (FillPatchIterator fpi(*oldht,Ydot,Ydot.nGrow(),cur_time,RhoYdot_Type,0,nspecies);
         fpi.isValid();
         ++fpi)
    {
        Ydot[fpi.index()].copy(fpi());
    }

    RhoH_to_Temp(get_new_data(State_Type));

    MultiFab& FuncCount = get_new_data(FuncCount_Type);

    for (FillPatchIterator fpi(*oldht,FuncCount,FuncCount.nGrow(),cur_time,FuncCount_Type,0,1);
         fpi.isValid();
         ++fpi)
    {
        FuncCount[fpi.index()].copy(fpi());
    }
}

//
// Inits the data on a new level that did not exist before regridding.
//

void
HeatTransfer::init ()
{
    NavierStokes::init();
 
    HeatTransfer& old      = getLevel(level-1);
    const Real    cur_time = old.state[State_Type].curTime();
    //
    // Get best ydot data.
    //
    FillCoarsePatch(get_new_data(RhoYdot_Type),0,cur_time,RhoYdot_Type,0,nspecies);

    RhoH_to_Temp(get_new_data(State_Type));

    if (new_T_threshold>0)
    {
        MultiFab& crse = old.get_new_data(State_Type);
        MultiFab& fine = get_new_data(State_Type);

        RhoH_to_Temp(crse,0); // Make sure T is current
        Real min_T_crse = crse.min(Temp);
        Real min_T_fine = min_T_crse * std::min(1.0, new_T_threshold);

        const int* ratio = crse_ratio.getVect();
        int Tcomp = (int)Temp;
        int Rcomp = (int)Density;
        int n_tmp = std::max( (int)Density, std::max( (int)Temp, std::max( last_spec, RhoH) ) );
        Array<Real> tmp(n_tmp);
        int num_cells_hacked = 0;
        for (MFIter mfi(fine); mfi.isValid(); ++mfi)
        {
            FArrayBox& fab = fine[mfi];
            const Box& box = mfi.validbox();
            num_cells_hacked += 
                FORT_CONSERVATIVE_T_FLOOR(box.loVect(), box.hiVect(),
                                          fab.dataPtr(), ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
                                          &min_T_fine, &Tcomp, &Rcomp, &first_spec, &last_spec, &RhoH,
                                          ratio, tmp.dataPtr(), &n_tmp);
        }

        ParallelDescriptor::ReduceIntSum(num_cells_hacked);

        if (num_cells_hacked > 0) {

            Real old_min = fine.min(Temp);
            RhoH_to_Temp(get_new_data(State_Type));
            Real new_min = fine.min(Temp);

            if (verbose && ParallelDescriptor::IOProcessor()) {
                std::cout << "...level data adjusted to reduce new extrema (" << num_cells_hacked
                          << " cells affected), new min = " << new_min << " (old min = " << old_min << ")" << '\n';
            }
        }
    }

    FillCoarsePatch(get_new_data(FuncCount_Type),0,cur_time,FuncCount_Type,0,1);
}

void
HeatTransfer::post_timestep (int crse_iteration)
{
    NavierStokes::post_timestep(crse_iteration);

    if (plot_reactions && level == 0)
    {
        const int Ndiag = auxDiag["REACTIONS"]->nComp();

        //
        // Multiply by the inverse of the coarse timestep.
        //
        const Real factor = 1.0 / crse_dt;

        for (int i = parent->finestLevel(); i >= 0; --i)
            getLevel(i).auxDiag["REACTIONS"]->mult(factor);

        for (int i = parent->finestLevel(); i > 0; --i)
        {
            HeatTransfer& clev = getLevel(i-1);
            HeatTransfer& flev = getLevel(i);

            MultiFab& Ydot_crse = *(clev.auxDiag["REACTIONS"]);
            MultiFab& Ydot_fine = *(flev.auxDiag["REACTIONS"]);

            NavierStokes::avgDown(clev.boxArray(),flev.boxArray(),
                                  Ydot_crse,Ydot_fine,
                                  i-1,i,0,Ndiag,parent->refRatio(i-1));
        }
    }
}

void
HeatTransfer::post_restart ()
{
    NavierStokes::post_restart();

    Real dummy  = 0;
    int MyProc  = ParallelDescriptor::MyProc();
    int step    = parent->levelSteps(0);
    int restart = 1;

    if (do_active_control)
        FORT_ACTIVECONTROL(&dummy,&dummy,&crse_dt,&MyProc,&step,&restart);

#ifdef PARTICLES
    if (level == 0)
    {
        BL_ASSERT(HTPC == 0);

        HTPC = new HTParticleContainer(parent);
        //
        // Make sure to call RemoveParticles() on exit.
        //
        BoxLib::ExecOnFinalize(RemoveParticles);

        HTPC->SetVerbose(pverbose);
        //
        // We want to be able to add new particles on a restart.
        // As well as the ability to write the particles out to an ascii file.
        //
        if (!restart_from_nonparticle_chkfile)
        {
            HTPC->Restart(parent->theRestartFile(), the_ht_particle_file_name);
        }

        if (!particle_restart_file.empty())
        {
            HTPC->InitFromAsciiFile(particle_restart_file,0);
        }

        if (!particle_output_file.empty())
        {
            HTPC->WriteAsciiFile(particle_output_file);
        }
    }
#endif
}

void
HeatTransfer::postCoarseTimeStep (Real cumtime)
{
    //
    // postCoarseTimeStep() is only called by level 0.
    //
    BL_ASSERT(level == 0);
}

void
HeatTransfer::post_regrid (int lbase,
                           int new_finest)
{
    NavierStokes::post_regrid(lbase, new_finest);
    //
    // FIXME: This may be necessary regardless, unless the interpolation
    //        to fine from coarse data preserves rho=sum(rho.Y)
    //
    if (do_set_rho_to_species_sum)
    {
        const int nGrow = 0;
        if (parent->levelSteps(0)>0 && level>lbase)
            set_rho_to_species_sum(get_new_data(State_Type),0,nGrow,0);
    }

#ifdef PARTICLES
    if (HTPC != 0)
    {
        HTPC->Redistribute();

        if (parent->finestLevel() > 0)
            HTPC->RemoveParticlesNotAtFinestLevel();
    }
#endif
}

void
HeatTransfer::checkPoint (const std::string& dir,
                          std::ostream&      os,
                          VisMF::How         how,
                          bool               dump_old)
{
    NavierStokes::checkPoint(dir,os,how,dump_old);

    if (level == 0)
    {
        if (ParallelDescriptor::IOProcessor())
        {
            const std::string tvfile = dir + "/" + typical_values_filename;
            std::ofstream tvos;
            tvos.open(tvfile.c_str(),std::ios::out);
            if (!tvos.good())
                BoxLib::FileOpenFailed(tvfile);
            Box tvbox(IntVect(),(NUM_STATE-1)*BoxLib::BASISV(0));
            int nComp = typical_values.size();
            FArrayBox tvfab(tvbox,nComp);
            for (int i=0; i<nComp; ++i) {
                tvfab.dataPtr()[i] = typical_values[i];
            }
            tvfab.writeOn(tvos);
        }

#ifdef PARTICLES
        if (HTPC != 0)
            HTPC->Checkpoint(dir,the_ht_particle_file_name);
#endif
    }
}

void
HeatTransfer::post_init (Real stop_time)
{
    if (level > 0)
        //
        // Nothing to sync up at level > 0.
        //
        return;

    const Real cur_time     = state[State_Type].curTime();
    const int  finest_level = parent->finestLevel();
    Real        dt_init     = 0.0;

    Array<Real> dt_save(finest_level+1);
    Array<int>  nc_save(finest_level+1);
    Real        dt_init2 = 0.0;
    Array<Real> dt_save2(finest_level+1);
    Array<int>  nc_save2(finest_level+1);
    //
    // Ensure state is consistent, i.e. velocity field satisfies initial
    // estimate of constraint, coarse levels are fine level averages, pressure
    // is zero.
    //
    post_init_state();
    //
    // Load typical values for each state component
    //
    set_typical_values(false);
    //
    // Estimate the initial timestepping.
    //
    post_init_estDT(dt_init,nc_save,dt_save,stop_time);
    //
    // Better estimate needs dt to estimate divu
    //
    const bool do_iter        = do_init_proj && projector;
    const int  init_divu_iter = do_iter ? num_divu_iters : 0;

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "doing num_divu_iters = " << num_divu_iters << '\n';

    for (int iter = 0; iter < init_divu_iter; ++iter)
    {
        //
        // Update species destruction rates in each level but not state.
        //
        if (nspecies > 0)
        {
            for (int k = 0; k <= finest_level; k++)
            {
                MultiFab& S_new = getLevel(k).get_new_data(State_Type);
                //
                // Don't update S_new in this strang_chem() call ...
                //
                MultiFab S_tmp(S_new.boxArray(),S_new.nComp(),0);

		MultiFab Forcing_tmp(S_new.boxArray(),nspecies+1,0);
		Forcing_tmp.setVal(0.);

                S_tmp.copy(S_new);  // Parallel copy

		getLevel(k).advance_chemistry(S_new,S_tmp,dt_save[k],Forcing_tmp,0);
            }
        }
        //
        // Recompute the velocity to obey constraint with chemistry and
        // divqrad and then average that down.
        //
        if (nspecies > 0)
        {
            for (int k = 0; k <= finest_level; k++)
            {
                MultiFab&  Divu_new = getLevel(k).get_new_data(Divu_Type);
                getLevel(k).calc_divu(cur_time,dt_save[k],Divu_new);
            }
            if (!hack_noavgdivu)
            {
                for (int k = finest_level-1; k >= 0; k--)
                {
                    HeatTransfer&   fine_lev = getLevel(k+1);
                    const BoxArray& fgrids   = fine_lev.grids;
                    
                    HeatTransfer&   crse_lev = getLevel(k);
                    const BoxArray& cgrids   = crse_lev.grids;
                    const IntVect&  fratio   = crse_lev.fine_ratio;

                    MultiFab& Divu_crse = crse_lev.get_new_data(Divu_Type);
                    MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);

                    crse_lev.NavierStokes::avgDown(cgrids,fgrids,
                                                   Divu_crse,Divu_fine,
                                                   k,k+1,0,1,fratio);
                }
            }
            //
            // Recompute the initial velocity field based on this new constraint
            //
            const Real divu_time = state[Divu_Type].curTime();

            int havedivu = 1;

            projector->initialVelocityProject(0,divu_time,havedivu);
            //
            // Average down the new velocity
            //
            for (int k = finest_level-1; k >= 0; k--)
            {
                const BoxArray& fgrids  = getLevel(k+1).grids;
                const BoxArray& cgrids  = getLevel(k).grids;
                MultiFab&       S_fine  = getLevel(k+1).get_new_data(State_Type);
                MultiFab&       S_crse  = getLevel(k).get_new_data(State_Type);
                IntVect&        fratio  = getLevel(k).fine_ratio;
                
                getLevel(k).NavierStokes::avgDown(cgrids,fgrids,
                                                  S_crse,S_fine,
                                                  k,k+1,Xvel,BL_SPACEDIM,
                                                  fratio);
            }
        }
        //
        // Estimate the initial timestepping again, using new velocity
        // (necessary?) Need to pass space to save dt, nc, but these are
        // hacked, just pass something.
        //
	// reset Ncycle to nref...
        //
	parent->setNCycle(nc_save);
        post_init_estDT(dt_init2, nc_save2, dt_save2, stop_time);
	//
	// Compute dt_init,dt_save as the minimum of the values computed
	// in the calls to post_init_estDT
	// Then setTimeLevel and dt_level to these values.
	//
	dt_init = std::min(dt_init,dt_init2);
	Array<Real> dt_level(finest_level+1,dt_init);

	parent->setDtLevel(dt_level);
	for (int k = 0; k <= finest_level; k++)
        {
	    dt_save[k] = std::min(dt_save[k],dt_save2[k]);
	    getLevel(k).setTimeLevel(cur_time,dt_init,dt_init);
        }
    }
    //
    // Initialize the pressure by iterating the initial timestep.
    //
    post_init_press(dt_init, nc_save, dt_save);
    //
    // Compute the initial estimate of conservation.
    //
    if (sum_interval > 0)
        sum_integrated_quantities();

}

void
HeatTransfer::sum_integrated_quantities ()
{
    const int finest_level = parent->finestLevel();
    const Real time        = state[State_Type].curTime();

    Real mass = 0.0;
    for (int lev = 0; lev <= finest_level; lev++)
	mass += getLevel(lev).volWgtSum("density",time);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "TIME= " << time << " MASS= " << mass;

    if (getChemSolve().index(fuelName) >= 0)
    {
        Real fuelmass = 0.0;
        std::string fuel = "rho.Y(" + fuelName + ")";
        for (int lev = 0; lev <= finest_level; lev++)
            fuelmass += getLevel(lev).volWgtSum(fuel,time);

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << " FUELMASS= " << fuelmass;

        int MyProc  = ParallelDescriptor::MyProc();
        int step    = parent->levelSteps(0);
        int restart = 0;

        if (do_active_control)
            FORT_ACTIVECONTROL(&fuelmass,&time,&crse_dt,&MyProc,&step,&restart);
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << '\n';

    Real rho_h    = 0.0;
    Real rho_temp = 0.0;

    for (int lev = 0; lev <= finest_level; lev++)
    {
        rho_temp += getLevel(lev).volWgtSum("rho_temp",time);
        rho_h     = getLevel(lev).volWgtSum("rhoh",time);
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "TIME= " << time << " LEV= " << lev << " RHOH= " << rho_h << '\n';
    }
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "TIME= " << time << " RHO*TEMP= " << rho_temp << '\n';

    if (verbose) {

        int old_prec = std::cout.precision(15);
        int nc=BL_SPACEDIM+1;
        Array<Real> vmin(nc), vmax(nc);
        Array<int> comps(nc);
        for (int n=0; n<BL_SPACEDIM; ++n) {
            comps[n] = Xvel+n;
        }
        comps[nc-1] = Temp;
        
        for (int lev = 0; lev <= finest_level; lev++)
        {
            MultiFab& newstate = getLevel(lev).get_data(State_Type,time);
            for (int n=0; n<=BL_SPACEDIM; ++n) {
                Real thismin = newstate.min(comps[n],0);
                Real thismax = newstate.max(comps[n],0);

                if (lev==0) {
                    vmin[n] = thismin;
                    vmax[n] = thismax;
                } else {
                    vmin[n] = std::min(thismin,vmin[n]);
                    vmax[n] = std::max(thismax,vmax[n]);
                }
            }
        }

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "min,max temp = " << vmin[nc-1] << ", " << vmax[nc-1] << '\n';
            for (int n=0; n<BL_SPACEDIM; ++n) {
                std::string str = (n==0 ? "xvel" : (n==1 ? "yvel" : "zvel") );
                std::cout << "min,max "<< str << "  = " << vmin[Xvel+n] << ", " << vmax[Xvel+n] << '\n';
            }
        }

        if (nspecies > 0)
        {
            Real min_sum, max_sum;
            for (int lev = 0; lev <= finest_level; lev++) {
                MultiFab* mf = getLevel(lev).derive("rhominsumrhoY",time,0);
                Real this_min = mf->min(0);
                Real this_max = mf->max(0);
                if (lev==0) {
                    min_sum = this_min;
                    max_sum = this_max;
                } else {
                    min_sum = std::min(this_min,min_sum);
                    max_sum = std::max(this_max,max_sum);
                }
                delete mf;
            }
            
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << "min,max rho-sum rho Y_l = "
                          << min_sum << ", " << max_sum << '\n';
            }
            
            for (int lev = 0; lev <= finest_level; lev++)
            {
                MultiFab* mf = getLevel(lev).derive("sumRhoYdot",time,0);
                Real this_min = mf->min(0);
                Real this_max = mf->max(0);
                if (lev==0) {
                    min_sum = this_min;
                    max_sum = this_max;
                } else {
                    min_sum = std::min(this_min,min_sum);
                    max_sum = std::max(this_max,max_sum);
                }
                delete mf;
            }
            
            if (ParallelDescriptor::IOProcessor()) {
                std::cout << "min,max sum RhoYdot = "
                          << min_sum << ", " << max_sum << '\n';
            }
        }       

        std::cout << std::setprecision(old_prec);
    }
}
    
void
HeatTransfer::post_init_press (Real&        dt_init,
			       Array<int>&  nc_save,
			       Array<Real>& dt_save)
{
    const int  nState          = desc_lst[State_Type].nComp();
    const int  nGrow           = 0;
    const Real cur_time        = state[State_Type].curTime();
    const int  finest_level    = parent->finestLevel();
    NavierStokes::initial_iter = true;
    //
    // Make space to save a copy of the initial State_Type state data
    //
    PArray<MultiFab> saved_state(finest_level+1,PArrayManage);

    if (init_iter > 0)
        for (int k = 0; k <= finest_level; k++)
            saved_state.set(k,new MultiFab(getLevel(k).grids,nState,nGrow));
    //
    // Iterate over the advance function.
    //
    for (int iter = 0; iter < init_iter; iter++)
    {
	//
	// Squirrel away copy of pre-advance State_Type state
	//
        for (int k = 0; k <= finest_level; k++)
	    MultiFab::Copy(saved_state[k],
                           getLevel(k).get_new_data(State_Type),
			   0,
                           0,
                           nState,
                           nGrow);

        for (int k = 0; k <= finest_level; k++ )
        {
            getLevel(k).advance(cur_time,dt_init,1,1);
        }
        //
        // This constructs a guess at P, also sets p_old == p_new.
        //
        MultiFab** sig = new MultiFab*[finest_level+1];

        for (int k = 0; k <= finest_level; k++)
        {
            sig[k] = getLevel(k).get_rho_half_time();
        }

        if (projector)
        {
            int havedivu = 1;
            projector->initialSyncProject(0,sig,parent->dtLevel(0),cur_time,
                                          havedivu);
        }
        delete [] sig;

        for (int k = finest_level-1; k>= 0; k--)
        {
            getLevel(k).avgDown();
        }
        for (int k = 0; k <= finest_level; k++)
        {
            //
            // Reset state variables to initial time, but
            // do not pressure variable, only pressure time.
            //
            getLevel(k).resetState(cur_time, dt_init, dt_init);
        }
	//
	// For State_Type state, restore state we saved above
	//
        for (int k = 0; k <= finest_level; k++)
	    MultiFab::Copy(getLevel(k).get_new_data(State_Type),
                           saved_state[k],
			   0,
                           0,
                           nState,
                           nGrow);

        NavierStokes::initial_iter = false;
    }

    if (init_iter <= 0)
        NavierStokes::initial_iter = false; // Just being compulsive -- rbp.

    NavierStokes::initial_step = false;
    //
    // Re-instate timestep.
    //
    for (int k = 0; k <= finest_level; k++)
    {
        getLevel(k).setTimeLevel(cur_time,dt_save[k],dt_save[k]);
        if (getLevel(k).state[Press_Type].descriptor()->timeType() ==
            StateDescriptor::Point)
          getLevel(k).state[Press_Type].setNewTimeLevel(cur_time+.5*dt_init);
    }
    parent->setDtLevel(dt_save);
    parent->setNCycle(nc_save);
}

//
// Reset the time levels to time (time) and timestep dt.
// This is done at the start of the timestep in the pressure
// iteration section.
//

void
HeatTransfer::resetState (Real time,
                          Real dt_old,
                          Real dt_new)
{
    NavierStokes::resetState(time,dt_old,dt_new);

    state[RhoYdot_Type].reset();
    state[RhoYdot_Type].setTimeLevel(time,dt_old,dt_new);

    state[FuncCount_Type].reset();
    state[FuncCount_Type].setTimeLevel(time,dt_old,dt_new);
}

void
HeatTransfer::avgDown ()
{
    if (level == parent->finestLevel()) return;

    HeatTransfer&   fine_lev = getLevel(level+1);
    const BoxArray& fgrids   = fine_lev.grids;
    MultiFab&       S_crse   = get_new_data(State_Type);
    MultiFab&       S_fine   = fine_lev.get_new_data(State_Type);

    NavierStokes::avgDown(grids,fgrids,S_crse,S_fine,
                          level,level+1,0,S_crse.nComp(),fine_ratio);
    //
    // Fill rho_ctime at the current and finer levels with the correct data.
    //
    for (int lev = level; lev <= parent->finestLevel(); lev++)
    {
        getLevel(lev).make_rho_curr_time();
    }
    //
    // Reset the temperature
    //
    RhoH_to_Temp(S_crse);
    //
    // Now average down pressure over time n-(n+1) interval.
    //
    MultiFab&       P_crse      = get_new_data(Press_Type);
    MultiFab&       P_fine_init = fine_lev.get_new_data(Press_Type);
    MultiFab&       P_fine_avg  = *fine_lev.p_avg;
    MultiFab&       P_fine      = initial_step ? P_fine_init : P_fine_avg;
    const BoxArray& P_fgrids    = fine_lev.state[Press_Type].boxArray();

    BoxArray crse_P_fine_BA = P_fgrids;

    crse_P_fine_BA.coarsen(fine_ratio);

    {
        MultiFab crse_P_fine(crse_P_fine_BA,1,0);

        for (MFIter mfi(P_fine); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            injectDown(crse_P_fine_BA[i],crse_P_fine[i],P_fine[i],fine_ratio);
        }

        P_crse.copy(crse_P_fine);  // Parallel copy
    }
    //
    // Next average down divu and dSdT at new time.
    //
    if (hack_noavgdivu) 
    {
        //
        // Now that state averaged down, recompute divu (don't avgDown,
        // since that will give a very different value, and screw up mac)
        //
        StateData& divuSD = get_state_data(Divu_Type);// should be const...
        const Real time   = divuSD.curTime();
        const Real dt     = time - divuSD.prevTime();
        calc_divu(time,dt,divuSD.newData());
    }
    else
    {
        MultiFab& Divu_crse = get_new_data(Divu_Type);
        MultiFab& Divu_fine = fine_lev.get_new_data(Divu_Type);
            
        NavierStokes::avgDown(grids,fgrids,Divu_crse,Divu_fine,
                              level,level+1,0,1,fine_ratio);
    }

    if (hack_noavgdivu)
    {
        //
        //  Now that have good divu, recompute dsdt (don't avgDown,
        //   since that will give a very different value, and screw up mac)
        //
        StateData& dsdtSD = get_state_data(Dsdt_Type);// should be const...
        MultiFab&  dsdt   = dsdtSD.newData();

        if (get_state_data(Divu_Type).hasOldData())
        {
            const Real time = dsdtSD.curTime();
            const Real dt   = time - dsdtSD.prevTime();
            calc_dsdt(time,dt,dsdt);
        }
        else
        {
            dsdt.setVal(0.0);
        }
    }
    else
    {
        MultiFab& Dsdt_crse = get_new_data(Dsdt_Type);
        MultiFab& Dsdt_fine = fine_lev.get_new_data(Dsdt_Type);
            
        NavierStokes::avgDown(grids,fgrids,Dsdt_crse,Dsdt_fine,
                              level,level+1,0,1,fine_ratio);
    }
}

void
HeatTransfer::differential_diffusion_update (MultiFab& Force,
                                             int       FComp,
                                             Real      theta,
                                             MultiFab& Dnew,
                                             int       DComp,
                                             MultiFab& DDnew,
                                             Real      theta_enthalpy)
{
    //
    // If theta_enthalpy > 0, recompute the D[Nspec+1] = Div(Fi.Hi) using
    // Fi.Hi based on solution here.  
    //
    BL_ASSERT(Force.boxArray() == grids);
    BL_ASSERT(FComp+Force.nComp()>=nspecies+1);
    BL_ASSERT(Dnew.boxArray() == grids);
    BL_ASSERT(DDnew.boxArray() == grids);
    BL_ASSERT(DComp+Dnew.nComp()>=nspecies+1);

    const Real strt_time = ParallelDescriptor::second();

    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! skipping spec diffusion " << '\n';

        for (int d = 0; d < BL_SPACEDIM; ++d)
        {
            SpecDiffusionFluxn[d]->setVal(0,0,nspecies+3);
            SpecDiffusionFluxnp1[d]->setVal(0,0,nspecies+3);
        }
        sumSpecFluxDotGradHn.setVal(0);
        sumSpecFluxDotGradHnp1.setVal(0);

        return;
    }
    
    MultiFab& S_old = get_old_data(State_Type);
    MultiFab& S_new = get_new_data(State_Type);
    Real prev_time = state[State_Type].prevTime();
    Real curr_time = state[State_Type].curTime();
    Real dt = curr_time - prev_time;
    BL_ASSERT(dt>0);

    const int nGrow = 0;

    if (theta == 0) // fully explicit update, assume all RHS terms in Force, zero np1 flux-based terms
    {
        MultiFab::Copy(S_new,Force,0,first_spec,nspecies+1,0);
        S_new.mult(dt,first_spec,nspecies);
        MultiFab::Add(S_new,S_old,first_spec,first_spec,nspecies+1,0);
        for (int d = 0; d < BL_SPACEDIM; ++d)
        {
            SpecDiffusionFluxnp1[d]->setVal(0,0,nspecies+3);
        }
        sumSpecFluxDotGradHnp1.setVal(0,0,1);
    }
    else
    {   
        MultiFab::Copy(S_new,S_old,first_spec,first_spec,nspecies+1,0);
        MultiFab* rho_half = 0;
        Diffusion::SolveMode solve_mode = Diffusion::ONEPASS;
        MultiFab **betanp1, **betan = 0; // Will not need betan since time-explicit pieces computed above

        diffusion->allocFluxBoxesLevel(betanp1,nGrow,nspecies+2); // fill rhoD, lambda/cp and lambda
        getDiffusivity(betanp1, curr_time, first_spec, 0, nspecies+1);
        getDiffusivity(betanp1, curr_time, Temp, nspecies+1, 1);
        //
        // Diffuse RhoY (and RhoH, unless theta_enthalpy>0 -- if so, will modify the source term for 
        //   RhoH with the results of the RhoY diffusion)
        //
        MultiFab* alpha = 0; // Never need alpha for RhoY, RhoH
        int alphaComp = 0;

        int nComp = (  theta_enthalpy > 0  ?  nspecies : nspecies+1 );
        for (int sigma = 0; sigma < nComp; ++sigma)
        {
            int betaComp = sigma;
            const int state_ind = first_spec + sigma;
            bool add_old_time_divFlux = false; // indicate that the rhs contains the time-explicit diff terms already
            int rho_flag = 2;
            diffusion->diffuse_scalar(dt,state_ind,theta,rho_half,rho_flag,
                                      SpecDiffusionFluxn,SpecDiffusionFluxnp1,sigma,&Force,sigma,alpha,
                                      alphaComp,betan,betanp1,betaComp,solve_mode,add_old_time_divFlux);

        }
        //
        // Remove theta scaling of fluxes, adjust so that the species fluxes sum to zero
        //  and compute - Div(Flux)
        //
        for (int d=0; d<BL_SPACEDIM; ++d) {
            (*SpecDiffusionFluxnp1[d]).mult(1/theta,0,nComp);
            (*SpecDiffusionFluxn[d]).mult(1/(1-theta),0,nComp);
        }
        //
        // Modify/update new-time fluxes to ensure sum of species fluxes = 0
        //
        bool grow_cells_already_filled = false;

        adjust_spec_diffusion_fluxes(curr_time,grow_cells_already_filled);

	// AJN FLUXREG
	// We have just performed the diffusion solve for Y_m.
	// If we are in the corrector, we have also performed the diffusion solve for h.
        // If updateFluxReg=T, we update VISCOUS flux registers:
	//   if we are in the predictor, ADD (1/2)*Gamma_{m,AD}^(0).
	//   if we are in the corrector, ADD Gamma_{m,AD}^(k+1),
	//                               ADD lambda^(k)/cp^(k) grad h_AD^(k+1).

	if (do_reflux && updateFluxReg)
	{
	  for (int d = 0; d < BL_SPACEDIM; d++)
	  {
	    for (MFIter fmfi(*SpecDiffusionFluxn[d]); fmfi.isValid(); ++fmfi)
	    {
	      if (level > 0)
	      {
		if (is_predictor)
		  getViscFluxReg().FineAdd((*SpecDiffusionFluxnp1[d])[fmfi],d,fmfi.index(),0,first_spec,nspecies  ,0.5*dt);
		else
		  getViscFluxReg().FineAdd((*SpecDiffusionFluxnp1[d])[fmfi],d,fmfi.index(),0,first_spec,nspecies+1,    dt);
	      }
	    }

	    if (level < parent->finestLevel())
	    {
	      if (is_predictor)
		getLevel(level+1).getViscFluxReg().CrseInit((*SpecDiffusionFluxnp1[d]),d,0,first_spec,nspecies  ,-0.5*dt,FluxRegister::ADD);
	      else
		getLevel(level+1).getViscFluxReg().CrseInit((*SpecDiffusionFluxnp1[d]),d,0,first_spec,nspecies+1,    -dt,FluxRegister::ADD);
	    }
	  }
	}

        flux_divergence(Dnew,DComp,SpecDiffusionFluxnp1,0,nComp,-1);

	if (theta_enthalpy > 0)
	{
            // update species with corrected fluxes
            // for species, snew = sold + dt*(theta*Dnp1 + F)
            for (MFIter mfi(Forcing); mfi.isValid(); ++mfi) 
            {
                const Box& box = mfi.validbox();
                const FArrayBox& f = Forcing[mfi];
                const FArrayBox& Sold = get_old_data(State_Type)[mfi];
                FArrayBox& Snew = get_new_data(State_Type)[mfi];
                const FArrayBox& Dnp1 = Dnew[mfi];
                
                Snew.copy(Dnp1,box,0,box,first_spec,nspecies);
                Snew.mult(theta,first_spec,nspecies);
                Snew.plus(f,box,box,0,first_spec,nspecies);
                Snew.mult(dt,first_spec,nspecies);
                Snew.plus(Sold,box,box,first_spec,first_spec,nspecies);
            }
            
            // compute enthalpy fluxes with correct species
            compute_enthalpy_fluxes(curr_time,betanp1,grow_cells_already_filled);

	    // AJN FLUXREG
	    // We have just computed the time-advanced "DD" terms for the forcing in the
	    //   predictor h implicit solve.
	    // If updateFluxReg=T, we update ADVECTIVE flux registers:
	    //   ADD (1/2)*h_m^n*(Gamma_{m,AD}^(0)-lambda^n/cp^n grad Y_{m,AD}^(0)).
	    if (do_reflux && updateFluxReg)
            {
	      for (int d = 0; d < BL_SPACEDIM; d++)
	      {
		for (MFIter fmfi(*SpecDiffusionFluxn[d]); fmfi.isValid(); ++fmfi)
		{
		  if (level > 0)
		  {
		    advflux_reg->FineAdd((*SpecDiffusionFluxn[d])[fmfi],d,fmfi.index(),nspecies+1,RhoH,1,0.5*dt);
		  }
		}

		if (level < parent->finestLevel())
		{
		  getAdvFluxReg(level+1).CrseInit((*SpecDiffusionFluxn[d]),d,nspecies+1,RhoH,1,-0.5*dt,FluxRegister::ADD);
		}
	      }
	    }

            flux_divergence(DDnew,0,SpecDiffusionFluxnp1,nspecies+1,1,-1);
            
            DDnew.mult(theta_enthalpy,0,1);
            MultiFab::Add(Force,DDnew,0,nspecies,1,0);
            DDnew.mult(1/theta_enthalpy,0,1);

            //
            // Now diffuse RhoH
            //
            int sigma = nspecies;
            int betaComp = sigma;
            const int state_ind = first_spec + sigma;
            bool add_old_time_divFlux = false; // indicate that the rhs contains the time-explicit diff terms already
            int rho_flag = 2;
            diffusion->diffuse_scalar(dt,state_ind,theta,rho_half,rho_flag,
                                      SpecDiffusionFluxn,SpecDiffusionFluxnp1,sigma,&Force,sigma,alpha,
                                      alphaComp,betan,betanp1,betaComp,solve_mode,add_old_time_divFlux);
            
            
            // Remove theta scaling of fluxes and compute - Div(lambda/cp Grad(H) )
            for (int d=0; d<BL_SPACEDIM; ++d) {
                (*SpecDiffusionFluxnp1[d]).mult(1/theta,sigma,1);
                (*SpecDiffusionFluxn[d]).mult(1/(1-theta),sigma,1);
            }

	    // AJN FLUXREG
	    // We have just finished the implicit solve for h in the predictor.
	    // If updateFluxReg=T, we update VISCOUS flux registers:
	    //  ADD (1/2)*lambda^n/cp^ngrad h_AD^(0)
	    if (do_reflux && updateFluxReg)
	    {
	      for (int d = 0; d < BL_SPACEDIM; d++)
		{
		  for (MFIter fmfi(*SpecDiffusionFluxn[d]); fmfi.isValid(); ++fmfi)
		  {
		    if (level > 0)
		      getViscFluxReg().FineAdd((*SpecDiffusionFluxnp1[d])[fmfi],d,fmfi.index(),nspecies,RhoH,1,0.5*dt);
		  }

		  if (level < parent->finestLevel())
		  {
		    getLevel(level+1).getViscFluxReg().CrseInit((*SpecDiffusionFluxnp1[d]),d,nspecies,RhoH,1,-0.5*dt,FluxRegister::ADD);
		  }
		}
	    }

            flux_divergence(Dnew,nspecies,SpecDiffusionFluxnp1,nspecies,1,-1);
        }

        diffusion->removeFluxBoxesLevel(betanp1);
    }
    
    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::differential_diffusion_update(): time: " << run_time << '\n';
    }
}

void
HeatTransfer::adjust_spec_diffusion_fluxes (Real                   time,
                                            bool                   grow_cells_already_filled)
{
    //
    // In this function we explicitly adjust the species diffusion fluxes so that their sum
    // is zero.  The fluxes are class member data, either SpecDiffusionFluxn or 
    // SpecDiffusionFluxnp1, depending on time
    //
    const TimeLevel whichTime = which_time(State_Type,time);
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    
    MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;

    MultiFab& S = get_data(State_Type,time);
    if (!grow_cells_already_filled)
    {
        // Fill grow cells in the state for RhoY, Temp
        // For Dirichlet physical boundary grow cells, this data will live on the cell face, otherwise
        // it will live at cell centers
        int nGrow = 1;
        BL_ASSERT(S.nGrow()>=nGrow);
        FillPatchIterator Tfpi(*this,S,nGrow,time,State_Type,Temp,1);
        FillPatchIterator Yfpi(*this,S,nGrow,time,State_Type,first_spec,nspecies);
        for ( ; Yfpi.isValid() && Tfpi.isValid(); ++Yfpi, ++Tfpi)
        {
            const Box& vbox = Yfpi.validbox();
            FArrayBox& fab = S[Yfpi];
            BoxList gcells = BoxLib::boxDiff(Box(vbox).grow(nGrow),vbox);
            for (BoxList::const_iterator it = gcells.begin(), end = gcells.end(); it != end; ++it)
            {
                const Box& gbox = *it;
                fab.copy(Tfpi(),gbox,0,gbox,Temp,1);
                fab.copy(Yfpi(),gbox,0,gbox,first_spec,nspecies);
            }
        }
    }

    // Get boundary info for Y (assume all Ys have the same boundary type
    const BCRec& Ybc = get_desc_lst()[State_Type].getBC(first_spec);

    // 
    // The following REPAIR_FLUX routine modifies the fluxes of all the species
    // to ensure that they sum to zero.  It requires the RhoY on valid + 1 grow (and, at least
    // as of this writing actually used the values of RhoY on edges that it gets by arithmetic
    // averaging.
    const Box& domain = geom.Domain();
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const Box& box   = mfi.validbox();
        FArrayBox& Y = S[mfi];
        int sCompY=first_spec;

        for (int d =0; d < BL_SPACEDIM; ++d)
        {
            FArrayBox& f = (*flux[d])[mfi.index()];
            FORT_REPAIR_FLUX(box.loVect(), box.hiVect(), domain.loVect(), domain.hiVect(),
                             f.dataPtr(),      ARLIM(f.loVect()),ARLIM(f.hiVect()),
                             Y.dataPtr(sCompY),ARLIM(Y.loVect()),ARLIM(Y.hiVect()),
                             &d, Ybc.vect());
        }
    }

}

void
HeatTransfer::compute_enthalpy_fluxes (Real                   time,
				       const MultiFab* const* beta,
                                       bool                   grow_cells_already_filled)
{


    BL_ASSERT(beta && beta[0]->nComp() == nspecies+2);

    const TimeLevel whichTime = which_time(State_Type,time);
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    
    MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;

    MultiFab& S = get_data(State_Type,time);

    MultiFab& sumSpecFluxDotGradH = (whichTime == AmrOldTime) ? sumSpecFluxDotGradHn : sumSpecFluxDotGradHnp1;
    //
    //  Compute species enthalpy on the edges, and increment the heat flux with the 
    //  flux of enthalpy due to species diffusion.  While here, we also compute the Fi.Grad(Hi) term
    //  required for the temperature equation.  Both the fluxes and the Fi.Grad(Hi) terms are stored
    //  in the class data.
    //
    const Box& domain = geom.Domain();
    const BCRec& Tbc = get_desc_lst()[State_Type].getBC(Temp);
    FArrayBox area[BL_SPACEDIM];

    // fill ghost cells for rhoY and Temp
    if (!grow_cells_already_filled)
      {
	FillPatchIterator rYfpi(*this,S,1,time,State_Type,first_spec,nspecies);
	FillPatchIterator Tfpi(*this,S,1,time,State_Type,Temp,1);

	for( ; rYfpi.isValid() && Tfpi.isValid(); ++rYfpi,++Tfpi)
	  {
	    S[rYfpi].copy(rYfpi(),0,first_spec,nspecies);
	    S[rYfpi].copy(Tfpi(),0,Temp,1);
	  }
      }

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const int i    = mfi.index();
        const Box& box = mfi.validbox();
            
        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            geom.GetFaceArea(area[dir],grids,i,dir,0);
            
        FArrayBox& fh = sumSpecFluxDotGradH[mfi];
        int FComp = 0;
            
        FArrayBox& T = S[mfi];
        int TComp = Temp;
            
        FArrayBox& RhoY = S[mfi];
        int RhoYComp = first_spec;
            
        int dComp = 0;
            
        const FArrayBox& rDx = (*beta[0])[i];
        FArrayBox& fix = (*flux[0])[i];
            
        const FArrayBox& rDy = (*beta[1])[i];
        FArrayBox& fiy = (*flux[1])[i];
            
#if BL_SPACEDIM == 3        
        const FArrayBox& rDz = (*beta[2])[i];
        FArrayBox& fiz = (*flux[2])[i];
#endif
        const Real* dx = geom.CellSize();
            
        FORT_ENTH_DIFF_TERMS(box.loVect(), box.hiVect(), domain.loVect(), domain.hiVect(), dx,
                             T.dataPtr(TComp), ARLIM(T.loVect()),  ARLIM(T.hiVect()),
                             RhoY.dataPtr(RhoYComp), ARLIM(RhoY.loVect()),  ARLIM(RhoY.hiVect()),
                                 
                             rDx.dataPtr(dComp),ARLIM(rDx.loVect()),ARLIM(rDx.hiVect()),
                             fix.dataPtr(FComp),ARLIM(fix.loVect()),ARLIM(fix.hiVect()),
                             area[0].dataPtr(), ARLIM(area[0].loVect()),ARLIM(area[0].hiVect()),
                                 
                             rDy.dataPtr(dComp),ARLIM(rDy.loVect()),ARLIM(rDy.hiVect()),
                             fiy.dataPtr(FComp),ARLIM(fiy.loVect()),ARLIM(fiy.hiVect()),
                             area[1].dataPtr(), ARLIM(area[1].loVect()),ARLIM(area[1].hiVect()),
#if BL_SPACEDIM == 3
                             rDz.dataPtr(dComp),ARLIM(rDz.loVect()),ARLIM(rDz.hiVect()),
                             fiz.dataPtr(FComp),ARLIM(fiz.loVect()),ARLIM(fiz.hiVect()),
                             area[2].dataPtr(), ARLIM(area[2].loVect()),ARLIM(area[2].hiVect()),
#endif
                             fh.dataPtr(),     ARLIM(fh.loVect()), ARLIM(fh.hiVect()),
                             Tbc.vect() );
    }
}
    
void
HeatTransfer::compute_OT_radloss (Real      time,
                                  int       nGrow,
                                  MultiFab& dqrad)
{
    //
    // Get optically thin radiation losses (+ve when energy LOST).
    //
    BL_ASSERT(do_OT_radiation || do_heat_sink);
    BL_ASSERT(nGrow <= dqrad.nGrow());
    BL_ASSERT(dqrad.boxArray() == grids);

    const Real* dx = geom.CellSize();

    Real p_amb, dpdt_factor;
    FORT_GETPAMB(&p_amb, &dpdt_factor);

    const Real Patm = p_amb / P1atm_MKS;
    const Real T_bg = 298.0;

    FArrayBox tmp;

    for (FillPatchIterator fpi(*this,dqrad,nGrow,time,State_Type,0,NUM_STATE);
         fpi.isValid();
         ++fpi)
    {
        FArrayBox& S   = fpi();
        const Box& box = S.box();

        tmp.resize(box,1);
        tmp.copy(S,Density,0,1);
        tmp.invert(1);
        for (int spec = first_spec; spec <= last_spec; spec++)
            S.mult(tmp,0,spec,1);

        FArrayBox& dqrad_fab = dqrad[fpi];

        tmp.resize(box,nspecies);

        getChemSolve().massFracToMoleFrac(tmp,S,box,first_spec,0);

        if (do_OT_radiation)
        {
            getChemSolve().getOTradLoss_TDF(dqrad_fab,S,tmp,Patm,T_bg,box,0,Temp,0);
        }
        else
        {
            dqrad_fab.setVal(0.0);
        }
        //
        // Add arbitrary heat sink.
        //
        if (do_heat_sink)
        {
            FArrayBox rad(box,1);
            FORT_RADLOSS(box.loVect(),box.hiVect(),
                         rad.dataPtr(),         ARLIM(rad.loVect()),ARLIM(rad.hiVect()),
                         S.dataPtr(Temp),       ARLIM(S.loVect()),  ARLIM(S.hiVect()),
                         S.dataPtr(first_spec), ARLIM(S.loVect()),  ARLIM(S.hiVect()),
                         dx, &Patm, &time);
            dqrad_fab.plus(rad,0,0,1);
        }
    }
}

void
HeatTransfer::diffuse_cleanup (MultiFab*&  delta_rhs,
                               MultiFab**& betan,
                               MultiFab**& betanp1,
                               MultiFab*&  alpha)
{
    delete delta_rhs;
    delete alpha;
    alpha = delta_rhs = 0;

    diffusion->removeFluxBoxesLevel(betan);
    diffusion->removeFluxBoxesLevel(betanp1);
}

void
HeatTransfer::diffuse_cleanup (MultiFab*&  delta_rhs,
                               MultiFab**& betan,
                               MultiFab**& betanp1)
{
    MultiFab* alpha = 0;
    diffuse_cleanup(delta_rhs, betan, betanp1, alpha);
}

void
HeatTransfer::velocity_diffusion_update (Real dt)
{
    //
    // Do implicit c-n solve for velocity
    // compute the viscous forcing
    // do following except at initial iteration--rbp, per jbb
    //
    if (is_diffusive[Xvel])
    {
        int rho_flag;
        if (do_mom_diff == 0)
        {
           rho_flag = 1;
        }
        else
        {
           rho_flag = 3;
        }

        MultiFab *delta_rhs = 0, **betan = 0, **betanp1 = 0;

        diffuse_velocity_setup(dt, delta_rhs, betan, betanp1);

        int rhsComp = 0;
        int betaComp = 0;
        diffusion->diffuse_velocity(dt,be_cn_theta,get_rho_half_time(),rho_flag,
                                    delta_rhs,rhsComp,betan,betanp1,betaComp);

        diffuse_cleanup(delta_rhs, betan, betanp1);
    }
}
    
void
HeatTransfer::diffuse_velocity_setup (Real        dt,
                                      MultiFab*&  delta_rhs,
                                      MultiFab**& betan,
                                      MultiFab**& betanp1)
{
    //
    // Do setup for implicit c-n solve for velocity.
    //
    BL_ASSERT(delta_rhs==0);
    BL_ASSERT(betan==0);
    BL_ASSERT(betanp1==0);
    const Real time = state[State_Type].prevTime();
    //
    // Assume always variable viscosity.
    //
    diffusion->allocFluxBoxesLevel(betan);
    diffusion->allocFluxBoxesLevel(betanp1);
    
    getViscosity(betan, time);
    getViscosity(betanp1, time+dt);
    //
    // Do the following only if it is a non-trivial operation.
    //
    if (S_in_vel_diffusion)
    {
        //
        // Include div mu S*I terms in rhs
        //  (i.e. make nonzero delta_rhs to add into RHS):
        //
        // The scalar and tensor solvers incorporate the relevant pieces of
        //  of Div(tau), provided the flow is divergence-free.  Howeever, if
        //  Div(U) =/= 0, there is an additional piece not accounted for,
        //  which is of the form A.Div(U).  For constant viscosity, Div(tau)_i
        //  = Lapacian(U_i) + mu/3 d[Div(U)]/dx_i.  For mu not constant,
        //  Div(tau)_i = d[ mu(du_i/dx_j + du_j/dx_i) ]/dx_i - 2mu/3 d[Div(U)]/dx_i
        //
        // As a convenience, we treat this additional term as a "source" in
        // the diffusive solve, computing Div(U) in the "normal" way we
        // always do--via a call to calc_divu.  This routine computes delta_rhs
        // if necessary, and stores it as an auxilliary rhs to the viscous solves.
        // This is a little strange, but probably not bad.
        //
        delta_rhs = new MultiFab(grids,BL_SPACEDIM,0);
        delta_rhs->setVal(0);
      
        MultiFab divmusi(grids,BL_SPACEDIM,0);
        //
	// Assume always variable viscosity.
        //
	diffusion->compute_divmusi(time,betan,divmusi);
	divmusi.mult((-2./3.)*(1.0-be_cn_theta),0,BL_SPACEDIM,0);
        (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
        //
	// Assume always variable viscosity.
        //
	diffusion->compute_divmusi(time+dt,betanp1,divmusi);
	divmusi.mult((-2./3.)*be_cn_theta,0,BL_SPACEDIM,0);
        (*delta_rhs).plus(divmusi,0,BL_SPACEDIM,0);
    }
}

void
HeatTransfer::getViscTerms (MultiFab& visc_terms,
                            int       src_comp, 
                            int       num_comp,
                            Real      time)
{
    //
    // Load "viscous" terms, starting from component = 0.
    //
    // JFG: for species, this procedure returns the *negative* of the divergence of 
    // of the diffusive fluxes.  specifically, in the mixture averaged case, the
    // diffusive flux vector for species k is
    //
    //       j_k = - rho D_k,mix grad Y_k
    //
    // so the divergence of the flux, div dot j_k, has a negative in it.  instead 
    // this procedure returns - div dot j_k to remove the negative.
    //
    // note the fluxes used in the code are extensive, that is, scaled by the areas
    // of the cell edges.  the calculation of the divergence is the sum of un-divided
    // differences of the extensive fluxes, all divided by volume.  so the effect is 
    // to give the true divided difference approximation to the divergence of the
    // intensive flux.
    //
    const int  last_comp = src_comp + num_comp - 1;
    int        icomp     = src_comp; // This is the current related state comp.
    int        load_comp = 0;        // Comp for result of current calculation.
    MultiFab** vel_visc  = 0;        // Potentially reused, raise scope
    const int  nGrow     = visc_terms.nGrow();
    //
    // Get Div(tau) from the tensor operator, if velocity and have non-const viscosity
    //
    if (src_comp < BL_SPACEDIM)
    {
        if (src_comp != Xvel || num_comp < BL_SPACEDIM)
            BoxLib::Error("tensor v -> getViscTerms needs all v-components at once");

        diffusion->allocFluxBoxesLevel(vel_visc);
        getViscosity(vel_visc, time);

        showMF("velVT",*viscn_cc,"velVT_viscn_cc",level);
        for (int dir=0; dir<BL_SPACEDIM; ++dir) {
            showMF("velVT",*(vel_visc[dir]),BoxLib::Concatenate("velVT_viscn_",dir,1),level);
        }

        int viscComp = 0;
        diffusion->getTensorViscTerms(visc_terms,time,vel_visc,viscComp);
        showMF("velVT",visc_terms,"velVT_visc_terms_1",level);

        icomp = load_comp = BL_SPACEDIM;
    }
    //
    // For Le != 1, need to get visc terms in a block
    //
    if (!unity_Le)
    {
	const int non_spec_comps = std::min(num_comp,
                                            std::max(0,first_spec-src_comp) + std::max(0,last_comp-last_spec));
	const bool has_spec = src_comp <= last_spec && last_comp >= first_spec;
	BL_ASSERT( !has_spec || (num_comp>=nspecies+non_spec_comps));
	if (has_spec)
	{
            if (do_mcdd)
            {
                if (load_comp!=0) {
                    BoxLib::Abort("Need to work out memory issues in getViscTerms for mcdd");
                }
                compute_mcdd_visc_terms(visc_terms,time,nGrow,DDOp::DD_Temp);
            }
	    else 
	    {
                //const int sCompY = first_spec - src_comp + load_comp;
		//compute_differential_diffusion_terms(visc_terms,sCompY,time,dt);
                BoxLib::Abort("Fix ME");
                showMF("velVT",visc_terms,"velVT_visc_terms_2",level);
            }
        }
    }
    //
    // Now, do the rest.
    // Get Div(visc.Grad(state)) or Div(visc.Grad(state/rho)), as appropriate.
    //
    for ( ; icomp <= last_comp; load_comp++, icomp++)
    {
	const bool is_spec  = icomp >= first_spec && icomp <= last_spec;
	if ( !(is_spec && !unity_Le) )
	{
	    if (icomp == Temp)
	    {
                if (do_mcdd) 
		{
		    // Do nothing, because in this case, was done above
		}
		else
		{
		  BoxLib::Abort("Should not call getTempViscTerms anymore, now terms are tied with species diffusion");
		}
	    }
	    else if (icomp == RhoH)
	    {
		if (do_mcdd) {
                    BoxLib::Abort("do we really want to get RhoH VT when do_mcdd?");
                }
                else
		{
		  BoxLib::Abort("Should not call getTempViscTerms anymore");
		}
	    }
	    else
	    {
		const int  rho_flag = Diffusion::set_rho_flag(diffusionType[icomp]);
		MultiFab** beta     = 0;

		if (icomp == Density)
                {
                    visc_terms.setVal(0.0,load_comp,1);
                }
                else
                {
                    //
                    // Assume always variable viscosity / diffusivity.
                    //
                    diffusion->allocFluxBoxesLevel(beta);
                    getDiffusivity(beta, time, icomp, 0, 1);

                    int betaComp = 0;
                    diffusion->getViscTerms(visc_terms,icomp-load_comp,
                                            icomp,time,rho_flag,beta,betaComp);
                    
                    diffusion->removeFluxBoxesLevel(beta);
                }
	    }
	}
    }
    //
    // Add Div(u) term if desired, if this is velocity, and if Div(u) is nonzero
    // If const-visc, term is mu.Div(u)/3, else it's -Div(mu.Div(u).I)*2/3
    //
    if (src_comp < BL_SPACEDIM && S_in_vel_diffusion)
    {
        if (num_comp < BL_SPACEDIM)
            BoxLib::Error("getViscTerms() need all v-components at once");
    
        MultiFab divmusi(grids,BL_SPACEDIM,1);
#ifndef NDEBUG
        showMF("velVT",get_old_data(Divu_Type),"velVT_divu",level);
#endif
        //
	// Assume always using variable viscosity.
        //
	diffusion->compute_divmusi(time,vel_visc,divmusi); // pre-computed visc above
	divmusi.mult((-2./3.),0,BL_SPACEDIM,0);
        showMF("velVT",divmusi,"velVT_divmusi",level);
        visc_terms.plus(divmusi,Xvel,BL_SPACEDIM,0);
        showMF("velVT",visc_terms,"velVT_visc_terms_3",level);
    }
    //
    // Clean up your mess ...
    //
    if (vel_visc)
        diffusion->removeFluxBoxesLevel(vel_visc);
    //
    // Ensure consistent grow cells
    //    
    if (nGrow > 0)
    {
        const int N = visc_terms.IndexMap().size();

#ifdef BL_USE_OMP
#pragma omp parallel for
#endif
        for (int i = 0; i < N; i++)
        {
            const int  k   = visc_terms.IndexMap()[i];
            FArrayBox& vt  = visc_terms[k];
            const Box& box = visc_terms.box(k);
            FORT_VISCEXTRAP(vt.dataPtr(),ARLIM(vt.loVect()),ARLIM(vt.hiVect()),
                            box.loVect(),box.hiVect(),&num_comp);
        }
        visc_terms.FillBoundary(0,num_comp);
        //
        // Note: this is a special periodic fill in that we want to
        // preserve the extrapolated grow values when periodic --
        // usually we preserve only valid data.  The scheme relies on
        // the fact that there is good data in the "non-periodic" grow cells.
        // ("good" data produced via VISCEXTRAP above)
        //
        geom.FillPeriodicBoundary(visc_terms,0,num_comp,true);
    }
}

void
HeatTransfer::fill_mcdd_boundary_data(Real time)
{
    const int nGrow = LinOp_grow;
    MultiFab S(grids,nspecies+1,nGrow);

    // Fill temp mf with Y,T
    const int compT = nspecies;
    const int compY = 0;

    for (FillPatchIterator T_fpi(*this,S,nGrow,time,State_Type,Temp,1),
             Y_fpi(*this,S,nGrow,time,State_Type,first_spec,nspecies),
             Rho_fpi(*this,S,nGrow,time,State_Type,Density,1);
         T_fpi.isValid() && Y_fpi.isValid() && Rho_fpi.isValid();
          ++T_fpi, ++Y_fpi, ++Rho_fpi)
    {
        // HACK: Assumes all fill-patched boundary data is coincident
        FArrayBox&       Sfab = S[T_fpi];
        const FArrayBox& Yfab = Y_fpi();
        const FArrayBox& Tfab = T_fpi();
        const FArrayBox& Rfab = Rho_fpi();

        Sfab.copy(Yfab,0,compY,nspecies);
        Sfab.copy(Tfab,0,compT,1);

        for (int n=0; n<nspecies; ++n)
            Sfab.divide(Rfab,0,compY+n,1);
    }

    BndryRegister* cbr = 0;
    const int compCT = nspecies;
    const int compCY = 0;

    if (level != 0)
    {
        BoxArray cgrids = BoxArray(grids).coarsen(crse_ratio);
        cbr = new BndryRegister();
        cbr->setBoxes(cgrids);
        for (OrientationIter fi; fi; ++fi)
            cbr->define(fi(),IndexType::TheCellType(),0,1,2,nspecies+1);

        const int nGrowC = 1;
        MultiFab SC(cgrids,nspecies+1,nGrowC,Fab_allocate);

        FillPatchIterator TC_fpi(parent->getLevel(level-1),SC,nGrowC,time,
                                 State_Type,Temp,1);
        FillPatchIterator YC_fpi(parent->getLevel(level-1),SC,nGrowC,time,
                                 State_Type,first_spec,nspecies);
        FillPatchIterator RhoC_fpi(parent->getLevel(level-1),SC,nGrowC,time,
                                   State_Type,Density,1);
        for ( ; TC_fpi.isValid() && YC_fpi.isValid() && RhoC_fpi.isValid();
              ++TC_fpi, ++YC_fpi, ++RhoC_fpi)
        {

            FArrayBox&       SCfab = SC[TC_fpi];
            const FArrayBox& YCfab = YC_fpi();
            const FArrayBox& TCfab = TC_fpi();
            const FArrayBox& RCfab = RhoC_fpi();

            SCfab.copy(YCfab,0,compY,nspecies);
            SCfab.copy(TCfab,0,compT,1);

            for (int n=0; n<nspecies; ++n)
                SCfab.divide(RCfab,0,compY+n,1);
        }

        cbr->copyFrom(SC,nGrowC,compY,compCY,nspecies);
        cbr->copyFrom(SC,nGrowC,compT,compCT,1);

    }

    int max_order = InterpBndryData::maxOrderDEF();
    MCDDOp.setBoundaryData(S,compT,S,compY,cbr,compCT,cbr,compCY,
                           get_desc_lst()[State_Type].getBC(Temp),
                           get_desc_lst()[State_Type].getBC(first_spec),max_order);

    if (level != 0) delete cbr;
}

void dumpProfileFab(const FArrayBox& fab,
                    const std::string file)
{
    const Box& box = fab.box();
    int imid = (int)(0.5*(box.smallEnd()[0] + box.bigEnd()[0]));
    IntVect iv1 = box.smallEnd(); iv1[0] = imid;
    IntVect iv2 = box.bigEnd(); iv2[0] = imid;
    iv1[1] = 115; iv2[1]=125;
    Box pb(iv1,iv2);
    std::cout << "dumping over " << pb << std::endl;

    std::ofstream osf(file.c_str());
    for (IntVect iv = pb.smallEnd(); iv <= pb.bigEnd(); pb.next(iv))
    {
        osf << iv[1] << " " << (iv[1]+0.5)*3.5/256 << " ";
        for (int n=0; n<fab.nComp(); ++n) 
            osf << fab(iv,n) << " ";
        osf << std::endl;
    }
    BoxLib::Abort();
}

void dumpProfile(const MultiFab& mf,
                 const std::string file)
{
    const FArrayBox& fab = mf[0];
    dumpProfileFab(fab,file);
}

void
HeatTransfer::compute_mcdd_visc_terms(MultiFab&           vtermsYH,
                                      Real                time,
                                      int                 nGrow,
                                      DDOp::DD_ApForTorRH whichApp,
                                      PArray<MultiFab>*   flux)
{
    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! zeroing spec+enth diffusion terms " << '\n';
        vtermsYH.setVal(0,0,nspecies+1,nGrow);
        return;            
    }

    BL_ASSERT(vtermsYH.boxArray()==grids);
    BL_ASSERT(vtermsYH.nComp() >= 1+nspecies); // Assume vt for Y start at 0 and for T/H at nspecies
    
    // If we werent given space to write fluxes, make some locally
    PArray<MultiFab>* tFlux;
    PArray<MultiFab> fluxLocal(BL_SPACEDIM,PArrayManage);
    if (flux==0)
    {
        for (int dir=0; dir<BL_SPACEDIM; ++dir)
            fluxLocal.set(dir,new MultiFab(BoxArray(grids).surroundingNodes(dir),nspecies+1,0));
        tFlux = &fluxLocal;
    }
    else
    {
       tFlux = flux;
    }
 
    // Load boundary data in operator (using c-f data where required)
    fill_mcdd_boundary_data(time);

    const int sCompY = 0;
    const int sCompT = nspecies;
    const int nGrowOp = 1;
    MultiFab S(grids,nspecies+1,nGrowOp);
    const MultiFab& State = get_data(State_Type,time);
    MultiFab::Copy(S,State,first_spec,sCompY,nspecies,0);
    MultiFab::Copy(S,State,Temp,sCompT,1,0);
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        FArrayBox&           Sfab = S[mfi];
        const FArrayBox& StateFab = State[mfi];
        const Box& box = mfi.validbox();
        for (int i=0; i<nspecies; ++i) {
            Sfab.divide(StateFab,box,Density,sCompY+i,1);
        }
    }
    showMF("mcddVT",S,"mcddVT_S");

    MCDDOp.setCoefficients(S,sCompT,S,sCompY);
    showMCDD("mcddOpVT",MCDDOp,"VT_DDOp");
    MCDDOp.applyOp(vtermsYH,S,*tFlux,whichApp);
    showMF("mcddVT",vtermsYH,"mcddVT_LofS");
    //
    // Add heat source terms to RhoH/T component
    //
    add_heat_sources(vtermsYH,nspecies,time,nGrow,1.0);
    showMF("mcddVT",vtermsYH,"mcddVT_addHS");

    //
    // Divide whole mess by cpmix at this time
    //   (WARNING: is almost can't be consistent with 
    //      any second order scheme for advancing T)
    //
    if (whichApp==DDOp::DD_Temp) 
    {
        FArrayBox cp;
        for (MFIter mfi(S); mfi.isValid(); ++mfi)
        {
            const Box& box = mfi.validbox();
            cp.resize(box,1);
            int sCompCp = 0;
            getChemSolve().getCpmixGivenTY(cp,S[mfi],S[mfi],box,sCompT,sCompY,sCompCp);
            vtermsYH[mfi].divide(cp,0,nspecies,1);
        }
    }
    showMF("mcddVT",vtermsYH,"mcddVT_divideCP");

    //
    // Ensure consistent grow cells
    //    
    if (nGrow > 0)
    {
        const int nc = nspecies+1;
        for (MFIter mfi(vtermsYH); mfi.isValid(); ++mfi)
        {
            FArrayBox& vtYH = vtermsYH[mfi];
            const Box& box = mfi.validbox();
            FORT_VISCEXTRAP(vtYH.dataPtr(),ARLIM(vtYH.loVect()),ARLIM(vtYH.hiVect()),
                            box.loVect(),box.hiVect(),&nc);
        }
        vtermsYH.FillBoundary(0,nspecies+1);
        //
        // Note: this is a special periodic fill in that we want to
        // preserve the extrapolated grow values when periodic --
        // usually we preserve only valid data.  The scheme relies on
        // the fact that there is computable data in the "non-periodic"
        // grow cells (produced via VISCEXTRAP above)
        //
        geom.FillPeriodicBoundary(vtermsYH,0,nspecies+1,true);
    }
    showMF("mcddVT",vtermsYH,"mcddVT_grow");
}

Real MFnorm(const MultiFab& mf,
            int             sComp,
            int             nComp)
{
    const int p = 0; // max norm
    Real normtot = 0.0;
    for (MFIter mfi(mf); mfi.isValid(); ++mfi)
        normtot = std::max(normtot,mf[mfi].norm(mfi.validbox(),p,sComp,nComp));
    ParallelDescriptor::ReduceRealMax(normtot);
    return normtot;
}

// Set T from H or H from T, depending on whichApp
static void
sync_H_T(MultiFab&                  S,
         const ChemDriver&          cd,
         const DDOp::DD_ApForTorRH& whichApp,
         int                        compY,
         int                        compT,
         int                        compH,
         int                        verbose,
         int                        level,
         const std::string&         id)
{
    if (whichApp == DDOp::DD_RhoH)
    {
        int  max_iters_taken = 0;
        bool error_occurred  = false;
        for (MFIter mfi(S); mfi.isValid(); ++mfi)
        {
            int iters_taken = cd.getTGivenHY(S[mfi],S[mfi],S[mfi],
                                             mfi.validbox(),compH,compY,compT);
            error_occurred = error_occurred  ||  iters_taken < 0;
            max_iters_taken = std::max(max_iters_taken, iters_taken);        
        }    
        ParallelDescriptor::ReduceBoolOr(error_occurred);
        ParallelDescriptor::ReduceIntMax(max_iters_taken);
        
        if (error_occurred && ParallelDescriptor::IOProcessor()) {
            std::string msg = "RhoH_to_Temp: error in H->T"+id;
            BoxLib::Error(msg.c_str());
        }
        if (verbose>2 && ParallelDescriptor::IOProcessor()) {
            std::string str = "\tmcdd_V::RhoH->T "+id+": max iters = ";
            levWrite(std::cout,level,str) << max_iters_taken << '\n';
        }
    }
    else
    {
        if (verbose>2 && ParallelDescriptor::IOProcessor()) {
            std::string str = "\tmcdd_V:: recomputing RhoH from T "+id;
            levWrite(std::cout,level,str) << '\n';
        }
        for (MFIter mfi(S); mfi.isValid(); ++mfi) {
            cd.getHmixGivenTY(S[mfi],S[mfi],S[mfi],
                              mfi.validbox(),compT,compY,compH);
        }    
    }
}

void
HeatTransfer::mcdd_fas_cycle(MCDD_MGParams&      p,
                             MultiFab&           S,
                             const MultiFab&     Rhs,
                             const MultiFab&     rhoCpInv,
                             PArray<MultiFab>&   flux,
                             Real                time,
                             Real                dt,
                             DDOp::DD_ApForTorRH whichApp)
{
    // This routine implements a V-cycle to relax the mcdd system, A(phi)=Rhs, 
    // where A(phi) = (rho.phi) - theta*dt*L(phi).  It is actually written as
    // Res = Rhs - A(phi) = 0.
    //
    // A couple notes about the implementation here:
    //
    // S contains Y, T, H and rho, in that order.
    // Rhs has nspecies+1 components, for the rhoY and rhoH/Temp equations, respectively.
    // flux is ordered identical to Rhs
    //
    // FIXME: next lines perhaps messy, since they have to agree with the calling routines setup...
    int sCompY = 0;
    int sCompT = nspecies;
    int sCompH = sCompT + 1;
    int sCompR = sCompH + 1;
    BL_ASSERT(sCompR<S.nComp());

    int dCompH = nspecies;

    const BoxArray& mg_grids = S.boxArray();

    // Form Res = Rhs - (rho.phi)^np1 - theta*dt*L(phi)^np1 
    int nGrow = 0;
    int nGrowOp = 1;
    int nComp = nspecies+1;

    Array<Real> mcdd_relax_factor(nComp,mcdd_relax_factor_Y);
    mcdd_relax_factor[nspecies] = mcdd_relax_factor_T;

    MultiFab T1(mg_grids,nComp,nGrow); // Leave name generic since it is used both for Lphi and Res
    MultiFab alpha(mg_grids,nComp,nGrowOp);

    int for_T0_H1 = (whichApp==DDOp::DD_Temp ? 0 : 1);
    ChemDriver& cd = getChemSolve();
    bool is_relaxed_nu1 = false;
    int mg_level = p.mg_level;

    std::string vgrp("mcdd");
    showMF(vgrp,S,"S_init",mg_level);
    showMF(vgrp,Rhs,"Rhs",mg_level);
    if (whichApp==DDOp::DD_Temp)
        showMF(vgrp,rhoCpInv,"rhoCpInv",mg_level);

    for (int iter=0; (iter<=p.nu1) && (p.status==HT_InProgress) && !is_relaxed_nu1; ++iter)
    {
        // Make T consistent with Y,H  FIXME: Seems to make sense only at level=0
        if (0 && mg_level==0 && p.iter==0 && iter==0) {
            sync_H_T(S,cd,whichApp,sCompY,sCompT,sCompH,mcdd_verbose,mg_level," (pre) ");
        }

        bool getAlpha = true;
        MCDDOp.applyOp(T1,S,flux,whichApp,mg_level,getAlpha,&alpha);
        if (whichApp==DDOp::DD_Temp) {
            MultiFab::Multiply(T1,rhoCpInv,0,dCompH,1,0);
            MultiFab::Multiply(alpha,rhoCpInv,0,dCompH,1,0);
        }
        showMF(vgrp,T1,"L_of_S",mg_level,iter);
        showMF(vgrp,alpha,"alpha",mg_level,iter);

        int res_only = (iter==mcdd_nu1[mg_level] ? 1 : 0);
        for (int i=0; i<nComp; ++i) {
            p.maxRes[i] = 0;
            p.maxCor[i] = 0;
        }

        for (MFIter mfi(S); mfi.isValid(); ++mfi)
        {
            FArrayBox& s = S[mfi];
            const FArrayBox& r = Rhs[mfi];
            const FArrayBox& L = T1[mfi];
            const FArrayBox& a = alpha[mfi];
            const Box& box = mfi.validbox();

            // Compute Res (returned as Lphi):
            //   Res = Rhs - A(S) = Rhs - rho.phi + theta.dt.Lphi 
            //       (or Res = Rhs - phi + theta.dt.Lphi for phi=T)
            // Then relax: phi = phi + f.Res/ALPHA, ALPHA=1+theta.dt.alpha
            //          since Lphi,alpha scaled by 1/rho.Cp already
            Real mult = -1;
            Real thetaDt = dt*be_cn_theta;
            FORT_HTDDRELAX(box.loVect(),box.hiVect(),
                           s.dataPtr(),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                           sCompY, sCompT, sCompH, sCompR,
                           L.dataPtr(),ARLIM(L.loVect()),ARLIM(L.hiVect()),
                           a.dataPtr(),ARLIM(a.loVect()),ARLIM(a.hiVect()),
                           r.dataPtr(),ARLIM(r.loVect()),ARLIM(r.hiVect()),
                           thetaDt, mcdd_relax_factor.dataPtr(), p.maxRes.dataPtr(), p.maxCor.dataPtr(),
                           for_T0_H1, res_only, mult);
        }
        ParallelDescriptor::ReduceRealMax(p.maxRes.dataPtr(),p.maxRes.size());
        ParallelDescriptor::ReduceRealMax(p.maxCor.dataPtr(),p.maxCor.size());
        showMF(vgrp,T1,"R_nu1",mg_level,iter);
        showMF(vgrp,S,"S_nu1",mg_level,iter);

        // Set (unscaled) initial residual
        if (p.iter==0 && iter==0) {
            for (int i=0; i<nComp; ++i) {
                p.maxRes_initial[i] = p.maxRes[i];
            }
#ifdef HT_SKIP_NITROGEN
            p.maxRes_initial[nspecies-1] = 1;
#endif
        }

        // Test if converged or stalled

        // Reduction of peak residual norm over all calls at this level
        int maxRes_redux_idx = 0;
        p.maxRes_norm = -1;
        Real peak_abs_res = -1.;
        int peak_abs_res_idx = 0;
        for (int i=0; i<nComp; ++i)
        {
            int tc = (i<nspecies ? first_spec+i : (whichApp==DDOp::DD_Temp ? Temp : RhoH) );
            Real thisRes = p.maxRes[i]/typical_values[tc];
            if (thisRes>peak_abs_res) {
                peak_abs_res = thisRes;
                peak_abs_res_idx = i;
            }
        }
        
        if (peak_abs_res<p.res_abs_tol)
        {
            p.status = HT_Solved;
            if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2)  {
                std::string str("mcdd - nu1: SOLVED - peak_abs_res (R1) < p.res_abs_tol (R2), I1 largest  R1,R2,I1:");
                levWrite(std::cout,mg_level,str)
                    << peak_abs_res << ", "
                    << p.res_abs_tol << ", "
                    << peak_abs_res_idx << std::endl;
            }
        }
        else
        {
            // Check residual reduction, but only for components with scaled mag of residual larger than abs tol
            for (int i=0; i<nComp; ++i)
            {
                int tc = (i<nspecies ? first_spec+i : (whichApp==DDOp::DD_Temp ? Temp : RhoH) );
                Real thisRes = p.maxRes[i]/typical_values[tc];
                if ( (thisRes>p.res_abs_tol)  && (p.maxRes_initial[i]/typical_values[tc]>p.res_abs_tol)) {
                    Real norm_res_i = std::abs(p.maxRes[i]/p.maxRes_initial[i]);
                    if (norm_res_i > p.maxRes_norm) {
                        p.maxRes_norm = norm_res_i;
                        maxRes_redux_idx = i;
                    }
                }
            }

            if (p.maxRes_norm<p.res_redux_tol) // Then all residuals bigger than abs_tol have been reduced by at least factor res_redux_tol
            {
                p.status = HT_Solved;
                if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2)  {
                    std::string str("mcdd - nu1: SOLVED - maxRes_norm(R1) < p.res_redux_tol (R2), I1 largest  R1,R2,I1:");
                    levWrite(std::cout,mg_level,str)
                        << p.maxRes_norm << ", "
                        << p.res_redux_tol << ", "
                        << maxRes_redux_idx << std::endl;
                }
            }
        }
        
#if 0   // NOT IMPLEMENTED SENSIBLY YET!!!

        // Reduction of peak residual norm over nu1 iterations
        if (p.status != HT_Solved)
        {
            if (iter==0) {
                iterRes_init = p.maxRes_norm;
            }
            else
            {
                if (p.maxRes_norm/iterRes_init < p.res_nu1_rtol) {
                    is_relaxed_nu1 = true;
                    
                    if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2)  {
                        std::string str("mcdd - nu1: SOLVED - p.maxRes_norm/iterRes_init (R1) < p.res_nu1_rtol (R2), I1 largest  R1,R2,I1:");
                        levWrite(std::cout,mg_level,str)
                            << p.maxRes_norm/iterRes_init << ", "
                            << p.res_nu1_rtol << ", "
                            << maxRes_redux_idx << std::endl;
                    }
                }
            }
        }
#endif
        
        // Magnitude of (scaled) correction during this relaxation sweep
        int max_cor_norm_idx = -1;
        p.maxCor_norm = 0;
        if (!res_only)
        {
            for (int i=0; i<nComp; ++i) {
                int tc = (i<nspecies ? first_spec+i : (whichApp==DDOp::DD_Temp ? Temp : RhoH) );
                if (p.maxCor[i]/typical_values[tc] > p.maxCor_norm) {
                    p.maxCor_norm = p.maxCor[i]/typical_values[tc];
                    max_cor_norm_idx = i;
                }
            }
            
            if (p.stalled_tol>0  && p.maxCor_norm<p.stalled_tol && p.status!=HT_Solved) {
                p.status = HT_Stalled;
                if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2)  {
                    std::string str("mcdd - nu1: STALLED - p.maxCor_norm (R1) < p.stalled_tol (R2)  R1,R2:");
                    levWrite(std::cout,mg_level,str)
                        << p.maxCor_norm << ", "
                        << p.stalled_tol << std::endl;
                }
            }
        }

        if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2) {
            if (mcdd_verbose < 6) 
            {
                std::string str("mcdd - nu1 (lev,iter,res,cor): (");
                levWrite(std::cout,mg_level,str)
                    << mg_level
                    << ", " << iter
                    << ", " << p.maxRes_norm;
                if (mcdd_verbose>3) {
                    std::cout << " [comp=" << maxRes_redux_idx;
                    if (mcdd_verbose>4) {
                        int tc = (peak_abs_res_idx<nspecies ? first_spec+peak_abs_res_idx : (whichApp==DDOp::DD_Temp ? Temp : RhoH) );
                        std::cout << ", Redux(R=" << p.maxRes[maxRes_redux_idx]
                                  << ", Rinit=" << p.maxRes_initial[maxRes_redux_idx]
                                  << ") Abs(comp=" << peak_abs_res_idx
                                  <<", R=" << p.maxRes[peak_abs_res_idx]/typical_values[tc]
                                  << ")";
                    }
                    std::cout << "]";
                }
                std::cout << ", " << p.maxCor_norm;
                if (mcdd_verbose>3) {
                    std::cout << " [comp=" << max_cor_norm_idx << "]";
                }
                std::cout << ")\n";
            }
            else
            {
                std::string str("mcdd - nu1 (lev,iter): (");
                levWrite(std::cout,mg_level,str)
                    << mg_level
                    << ", " << iter << ")" << endl;
                str = "                         ";
                levWrite(std::cout,mg_level,str)
                    << '\t' << "Residual   " 
                    << '\t' << "Res/typ    " 
                    << '\t' << "Res init   " 
                    << '\t' << "R_init/typ " 
                    << '\t' << "Correction "
                    << '\t' << "Corr/typ   "
                    << '\t' << "TypicalVal "
                    << endl; 
                for (int i=0; i<nComp; ++i) 
                {
                    int tc = (i<nspecies ? first_spec+i : (whichApp==DDOp::DD_Temp ? Temp : RhoH) );
                    levWrite(std::cout,mg_level,str) << i
                                                     << '\t' << p.maxRes[i]
                                                     << '\t' <<  p.maxRes[i]/typical_values[tc]
                                                     << '\t' <<  p.maxRes_initial[i]
                                                     << '\t' <<  p.maxRes_initial[i]/typical_values[tc]
                                                     << '\t' <<  p.maxCor[i]
                                                     << '\t' <<  p.maxCor[i]/typical_values[tc]
                                                     << '\t' <<  typical_values[tc]
                                                     << endl;
                }
                levWrite(std::cout,mg_level,str) << '\t' << "Index of maximum normalized residual: " << peak_abs_res_idx << endl; 
                levWrite(std::cout,mg_level,str) << '\t' << "Index of maximum residual reduction: " << maxRes_redux_idx << endl; 
                if (max_cor_norm_idx>=0)
                    levWrite(std::cout,mg_level,str) << '\t' << "Index of maximum normalized correction: " << max_cor_norm_idx << endl << endl; 
                else
                    levWrite(std::cout,mg_level,str) << '\t' << "residual eval only - no correction" << endl<< endl; 
            }
        }
    }  // nu1 iters

    if ((p.status==HT_InProgress) && p.num_coarser>0) {

        // FIXME: BA calcs + allocs can be done ahead of time
        const IntVect MGIV(D_DECL(2,2,2));
        const BoxArray c_grids = BoxArray(mg_grids).coarsen(MGIV);
        MultiFab SC(c_grids,nspecies+3,nGrowOp);
        MultiFab T1C(c_grids,nspecies+1,nGrow); // Will be used for LphiC, ResC
        MultiFab T2C(c_grids,nspecies+3,nGrow); // Will be used for RhsC and SC_save, need extra comps for Rho and H        
        MultiFab rhoCpInvC(c_grids,1,nGrow);
        PArray<MultiFab> fluxC(BL_SPACEDIM,PArrayManage);
        for (int dir=0; dir<BL_SPACEDIM; ++dir)
            fluxC.set(dir,new MultiFab(BoxArray(c_grids).surroundingNodes(dir),
                                       nspecies+1,0));

        //
        // Build Rhs for the coarse problem = Avg(Res) + Op(Avg(S))
        //  (here, construct Avg(Res), then use HTDDRELAX to do Op and add)
        //
        if (whichApp==DDOp::DD_Temp)
        {
            const_cast<MultiFab&>(rhoCpInv).invert(1.,0,1,0); // Will un-invert just below...not really changing it...much
            DDOp::average(rhoCpInvC,0,rhoCpInv,0,1);  // Avg(rho.cp)
            const_cast<MultiFab&>(rhoCpInv).invert(1.,0,1,0);
            rhoCpInvC.invert(1.,0,1,0);
            showMF(vgrp,rhoCpInvC,"Avg_RhoCpInvC",mg_level);
        }

        // HACK
        T2C.setVal(0,nspecies+1,2); // top two components arent used, but the NaNs irritate amrvis

        T2C.setVal(0);

        for (MFIter mfi(T1); mfi.isValid(); ++mfi)
        {
            BL_ASSERT(!T1[mfi].contains_nan(mfi.validbox(),0,nspecies+1));
        }

        DDOp::average(T2C,0,T1,0,nspecies+1);  // Avg(Res)=Avg(T1)=T2C

        for (MFIter mfi(T2C); mfi.isValid(); ++mfi)
        {
            BL_ASSERT(!T2C[mfi].contains_nan(mfi.validbox(),0,nspecies+1));
        }

        showMF(vgrp,T2C,"Avg_Res",mg_level);

        // Rho-weight the state, then average (leave fine data rho-weighted for now...)
        for (MFIter mfi(S); mfi.isValid(); ++mfi) {
            FArrayBox& Sfab = S[mfi];
            for (int n=0; n<nspecies; ++n) {
                Sfab.mult(Sfab,sCompR,sCompY+n,1);
            }
            Sfab.mult(Sfab,sCompR,sCompH,1);
        }

        for (MFIter mfi(S); mfi.isValid(); ++mfi)
        {
            BL_ASSERT(!S[mfi].contains_nan(mfi.validbox(),0,nspecies+3));
        }

        DDOp::average(SC,0,S,0,nspecies+3); // Includes generation of rho, rhoh on coarse

        for (MFIter mfi(SC); mfi.isValid(); ++mfi)
        {
            BL_ASSERT(!SC[mfi].contains_nan(mfi.validbox(),0,nspecies+3));
        }

        showMF(vgrp,SC,"Avg_Srho",mg_level);

        // Unweight coarse state
        for (MFIter mfi(SC); mfi.isValid(); ++mfi) {
            FArrayBox& SCfab = SC[mfi];
            SCfab.invert(1.,sCompR,1);
            for (int n=0; n<nspecies; ++n) {
                SCfab.mult(SCfab,sCompR,sCompY+n,1);
            }
            SCfab.mult(SCfab,sCompR,sCompH,1);
            SCfab.invert(1.,sCompR,1);
        }
        showMF(vgrp,SC,"Avg_S",mg_level);

        MCDD_MGParams pC(p);
        pC.mg_level++;
        pC.status = HT_InProgress;

        pC.num_coarser--;
        if (pC.num_coarser<1  &&  mcdd_nub>0) {
            pC.nu1 = mcdd_nub;
            pC.nu2 = 0;
        } else {
            pC.nu1 = mcdd_nu1[pC.mg_level];
            pC.nu2 = mcdd_nu2[pC.mg_level];
        }

        bool getAlphaC = false;
        MultiFab* alphaC = 0;
        MCDDOp.applyOp(T1C,SC,fluxC,whichApp,pC.mg_level,getAlphaC,alphaC);
        if (whichApp==DDOp::DD_Temp) {
            MultiFab::Multiply(T1,rhoCpInv,0,dCompH,1,0);
            MultiFab::Multiply(alpha,rhoCpInv,0,dCompH,1,0);
        }
        showMF(vgrp,T1C,"L_of_Avg_S",mg_level);

        int res_only = 1; // here, using HTDDRELAX to build (RhsC + Op(SC)) only
        for (MFIter mfi(SC); mfi.isValid(); ++mfi)
        {
            FArrayBox& sC = SC[mfi];
            const FArrayBox& rC = T2C[mfi]; // Averaged fine grid residual
            const FArrayBox& lC = T1C[mfi]; // Op(Avg(S)) = Op(SC)
            const FArrayBox& aC = alpha[mfi]; // unused, but pass something
            const Box& box = mfi.validbox();

            // Compute Rhs to the coarse problem = Avg(Rhs) + Op(S_avg) [result in T1C]
            Real mult = 1;
            Real thetaDt = be_cn_theta * dt;
            FORT_HTDDRELAX(box.loVect(),box.hiVect(),
                           sC.dataPtr(),ARLIM(sC.loVect()),ARLIM(sC.hiVect()),
                           sCompY, sCompT, sCompH, sCompR,
                           lC.dataPtr(),ARLIM(lC.loVect()),ARLIM(lC.hiVect()),
                           aC.dataPtr(),ARLIM(aC.loVect()),ARLIM(aC.hiVect()),
                           rC.dataPtr(),ARLIM(rC.loVect()),ARLIM(rC.hiVect()),
                           thetaDt, mcdd_relax_factor.dataPtr(), pC.maxRes.dataPtr(), pC.maxCor.dataPtr(),
                           for_T0_H1, res_only, mult);

        }
        ParallelDescriptor::ReduceRealMax(pC.maxRes.dataPtr(),pC.maxRes.size());
        ParallelDescriptor::ReduceRealMax(pC.maxCor.dataPtr(),pC.maxCor.size());

        showMF(vgrp,T1C,"Res_of_Avg_S",mg_level);

        // Save a copy of pre-relaxed coarse state, in T2C
        MultiFab::Copy(T2C,SC,0,0,nspecies+3,0);

        // Do coarser relax and ignore solver status
        for (pC.iter=0; pC.iter<p.gamma; ++pC.iter) {
            mcdd_fas_cycle(pC,SC,T1C,rhoCpInvC,fluxC,time,dt,whichApp);
        }

        showMF(vgrp,SC,"SC_postV",mg_level);

        // Compute increment
        MultiFab::Subtract(SC,T2C,sCompY,sCompY,nspecies+2,nGrow);
        showMF(vgrp,SC,"eC",mg_level);

        // Rho-weight the increment to the state
        for (MFIter mfi(SC); mfi.isValid(); ++mfi) {
            FArrayBox& SCfab = SC[mfi];
            for (int n=0; n<nspecies; ++n) {
                SCfab.mult(SCfab,sCompR,sCompY+n,1);
            }
            SCfab.mult(SCfab,sCompR,sCompH,1);
        }

        showMF(vgrp,SC,"eCrho",mg_level);

        // Increment (already rho-weighted) fine solution with rho-weighted
        //   interpolated correction, then unweight
        DDOp::interpolate(S,0,SC,0,nspecies+2); // S += I(dSC)

        showMF(vgrp,S,"Srho_plus_inc",mg_level);

        // Unweight the fine state
        for (MFIter mfi(S); mfi.isValid(); ++mfi) {
            FArrayBox& Sfab = S[mfi];
            Sfab.invert(1.,sCompR,1);
            for (int n=0; n<nspecies; ++n) {
                Sfab.mult(Sfab,sCompR,sCompY+n,1);
            }
            Sfab.mult(Sfab,sCompR,sCompH,1);
            Sfab.invert(1.,sCompR,1);
        }

        if (0 && mg_level==0) {
            sync_H_T(S,cd,whichApp,sCompY,sCompT,sCompH,mcdd_verbose,mg_level," (post coarse) ");
        }

        showMF(vgrp,S,"S_pre_nu2",mg_level);

        // Do post smooth (if not at the bottom)
        bool is_relaxed_nu2 = false;
        for (int iter=0; (iter<=p.nu2) && (p.status==HT_InProgress) && !is_relaxed_nu2; ++iter)
        {
            // Make T consistent with Y,H
            if (0 && mg_level==0 && iter==0) {
                sync_H_T(S,cd,whichApp,sCompY,sCompT,sCompH,mcdd_verbose,mg_level," (post) ");
            }
            
            // Compute Res = Rhs - A(S) = Rhs - rho.phi + theta.dt.Lphi
            // and then relax: phi = phi + f.Res/alpha
            bool getAlpha = true;
            MCDDOp.applyOp(T1,S,flux,whichApp,mg_level,getAlpha,&alpha);
            if (whichApp==DDOp::DD_Temp) {
                MultiFab::Multiply(T1,rhoCpInv,0,dCompH,1,0);
                MultiFab::Multiply(alpha,rhoCpInv,0,dCompH,1,0);
            }
            showMF(vgrp,T1,"L_of_S2",mg_level,iter);
            showMF(vgrp,alpha,"alpha2",mg_level,iter);

            for (int i=0; i<nComp; ++i) {
                p.maxCor[i] = 0;
                p.maxRes[i] = 0;
            }
            
            int res_only = (iter==mcdd_nu2[mg_level] ? 1 : 0);
            for (MFIter mfi(S); mfi.isValid(); ++mfi)
            {
                FArrayBox& s = S[mfi];
                const FArrayBox& r = Rhs[mfi];
                const FArrayBox& L = T1[mfi];
                const FArrayBox& a = alpha[mfi];
                const Box& box = mfi.validbox();
                
                // Compute Res (returned as Lphi):
                //   Res = Rhs - A(S) = Rhs - rho.phi + theta.dt.Lphi (or Rhs - phi + theta.dt.Lphi if phi=Temp)
                // Then relax: phi = phi + f.Res/alpha
                Real mult = -1;
                Real thetaDt = be_cn_theta * dt;
                FORT_HTDDRELAX(box.loVect(),box.hiVect(),
                               s.dataPtr(),ARLIM(s.loVect()),ARLIM(s.hiVect()),
                               sCompY, sCompT, sCompH, sCompR,
                               L.dataPtr(),ARLIM(L.loVect()),ARLIM(L.hiVect()),
                               a.dataPtr(),ARLIM(a.loVect()),ARLIM(a.hiVect()),
                               r.dataPtr(),ARLIM(r.loVect()),ARLIM(r.hiVect()),
                               thetaDt, mcdd_relax_factor.dataPtr(), p.maxRes.dataPtr(), p.maxCor.dataPtr(),
                               for_T0_H1, res_only, mult);
            } // S mfi
            
            ParallelDescriptor::ReduceRealMax(p.maxCor.dataPtr(),p.maxCor.size());
            ParallelDescriptor::ReduceRealMax(p.maxRes.dataPtr(),p.maxRes.size());

            showMF(vgrp,T1,"R_nu2",mg_level,iter);
            showMF(vgrp,S,"S_nu2",mg_level,iter);

            // Test if converged or stalled
            
            // Reduction of peak residual norm over all calls at this level
            int maxRes_redux_idx = 0;
            p.maxRes_norm = -1;
            Real peak_abs_res = -1.;
            int peak_abs_res_idx = 0;
            for (int i=0; i<nComp; ++i)
            {
                int tc = (i<nspecies ? first_spec+i : (whichApp==DDOp::DD_Temp ? Temp : RhoH) );
                Real thisRes = p.maxRes[i]/typical_values[tc];
                if (thisRes>peak_abs_res) {
                    peak_abs_res = thisRes;
                    peak_abs_res_idx = i;
                }
            }
            
            if (peak_abs_res<p.res_abs_tol)
            {
                p.status = HT_Solved;
                if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2)  {
                    std::string str("mcdd - nu2: SOLVED - peak_abs_res (R1) < p.res_abs_tol (R2), I1 largest  R1,R2,I1:");
                    levWrite(std::cout,mg_level,str)
                        << peak_abs_res << ", "
                        << p.res_abs_tol << ", "
                        << peak_abs_res_idx << std::endl;
                }
            }
            else
            {
                // Check residual reduction, but only for components with scaled mag of residual larger than abs tol
                for (int i=0; i<nComp; ++i)
                {
                    int tc = (i<nspecies ? first_spec+i : (whichApp==DDOp::DD_Temp ? Temp : RhoH) );
                    Real thisRes = p.maxRes[i]/typical_values[tc];
                    if ( (thisRes>p.res_abs_tol)  && (p.maxRes_initial[i]/typical_values[tc]>p.res_abs_tol)) {
                        Real norm_res_i = std::abs(p.maxRes[i]/p.maxRes_initial[i]);
                        if (norm_res_i > p.maxRes_norm) {
                            p.maxRes_norm = norm_res_i;
                            maxRes_redux_idx = i;
                        }
                    }
                }
                
                if (p.res_abs_tol>0  &&  (peak_abs_res<p.res_abs_tol)) 
                {
                    p.status = HT_Solved;
                    if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2)  {
                        std::string str("mcdd - nu2: SOLVED - peak_abs_res (R1) < p.res_abs_tol (R2), I1 largest  R1,R2,I1:");
                        levWrite(std::cout,mg_level,str)
                            << peak_abs_res << ", "
                            << p.res_abs_tol << ", "
                            << peak_abs_res_idx << std::endl;
                    }
                }
                else
                {
                    if (p.maxRes_norm<p.res_redux_tol) // Then have we reduced all by factor res_redux_tol
                    {
                        p.status = HT_Solved;
                        if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2)  {
                            std::string str("mcdd - nu2: SOLVED - maxRes_norm(R1) < p.res_redux_tol (R2), I1 largest  R1,R2,I1:");
                            levWrite(std::cout,mg_level,str)
                                << p.maxRes_norm << ", "
                                << p.res_redux_tol << ", "
                                << maxRes_redux_idx << std::endl;
                        }
                    }
                }
            }

#if 0   // NOT IMPLEMENTED SENSIBLY YET!!!

            // Reduction of peak residual norm over nu2 iterations
            if (p.status != HT_Solved)
            {
                if (iter==0) {
                    iterRes_init = p.maxRes_norm;
                }
                else
                {
                    if (p.maxRes_norm/iterRes_init < p.res_nu2_rtol) {
                        is_relaxed_nu2 = true;

                        if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2)  {
                            std::string str("mcdd - nu2: SOLVED - p.maxRes_norm/iterRes_init (R1) < p.res_nu2_rtol (R2), I1 largest  R1,R2,I1:");
                            levWrite(std::cout,mg_level,str)
                                << p.maxRes_norm/iterRes_init << ", "
                                << p.res_nu2_rtol << ", "
                                << maxRes_redux_idx << std::endl;
                        }
                    }
                }
            }
#endif

            // Magnitude of (scaled) correction during this relaxation sweep
            int max_cor_norm_idx = -1;
            p.maxCor_norm = 0;
            if (!res_only)
            {
                for (int i=0; i<nComp; ++i) {
                    int tc = (i<nspecies ? first_spec+i : (whichApp==DDOp::DD_Temp ? Temp : RhoH) );
                    if (p.maxCor[i]/typical_values[tc] > p.maxCor_norm) {
                        p.maxCor_norm = p.maxCor[i]/typical_values[tc];
                        max_cor_norm_idx = i;
                    }
                }
            
                if (p.stalled_tol>0  && p.maxCor_norm<p.stalled_tol && p.status!=HT_Solved) {
                    p.status = HT_Stalled;
                    if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2)  {
                        std::string str("mcdd - nu2: STALLED - p.maxCor_norm (R1) < p.stalled_tol (R2)  R1,R2:");
                        levWrite(std::cout,mg_level,str)
                            << p.maxCor_norm << ", "
                            << p.stalled_tol << std::endl;
                    }
                }
            }
                
            if (ParallelDescriptor::IOProcessor() && mcdd_verbose>2) {
                if (mcdd_verbose < 6) 
                {
                    std::string str("mcdd - nu2: (lev,iter,res,cor): (");
                    levWrite(std::cout,mg_level,str)
                        << mg_level
                        << ", " << iter
                        << ", " << p.maxRes_norm;
                    if (mcdd_verbose>3) {
                        std::cout << " [comp=" << maxRes_redux_idx;
                        if (mcdd_verbose>4) {
                            int tc = (peak_abs_res_idx<nspecies ? first_spec+peak_abs_res_idx : (whichApp==DDOp::DD_Temp ? Temp : RhoH) );
                            std::cout << ", Redux(R=" << p.maxRes[maxRes_redux_idx]
                                      << ", Rinit=" << p.maxRes_initial[maxRes_redux_idx]
                                      << ") Abs(comp=" << peak_abs_res_idx
                                      <<", R=" << p.maxRes[peak_abs_res_idx]/typical_values[tc]
                                      << ")";
                        }
                        std::cout << "]";
                    }
                    std::cout << ", " << p.maxCor_norm;
                    if (mcdd_verbose>3) {
                        std::cout << " [comp=" << max_cor_norm_idx << "]";
                    }
                    std::cout << ")\n";
                }
                else
                {
                    std::string str("mcdd - nu1 (lev,iter): (");
                    levWrite(std::cout,mg_level,str)
                        << mg_level
                        << ", " << iter << ")" << endl;
                    str = "                         ";
                    levWrite(std::cout,mg_level,str)
                        << '\t' << "Residual   " 
                        << '\t' << "Res/typ    " 
                        << '\t' << "Res init   " 
                        << '\t' << "R_init/typ " 
                        << '\t' << "Correction "
                        << '\t' << "Corr/typ   "
                        << '\t' << "TypicalVal "
                        << endl; 
                    for (int i=0; i<nComp; ++i) 
                    {
                        int tc = (i<nspecies ? first_spec+i : (whichApp==DDOp::DD_Temp ? Temp : RhoH) );
                        levWrite(std::cout,mg_level,str) << i
                                                         << '\t' << p.maxRes[i]
                                                         << '\t' <<  p.maxRes[i]/typical_values[tc]
                                                         << '\t' <<  p.maxRes_initial[i]
                                                         << '\t' <<  p.maxRes_initial[i]/typical_values[tc]
                                                         << '\t' <<  p.maxCor[i]
                                                         << '\t' <<  p.maxCor[i]/typical_values[tc]
                                                         << '\t' <<  typical_values[tc]
                                                         << endl;
                    }
                    levWrite(std::cout,mg_level,str) << '\t' << "Index of maximum normalized residual: " << peak_abs_res_idx << endl; 
                    levWrite(std::cout,mg_level,str) << '\t' << "Index of maximum residual reduction: " << maxRes_redux_idx << endl; 
                    if (max_cor_norm_idx>=0)
                        levWrite(std::cout,mg_level,str) << '\t' << "Index of maximum normalized correction: " << max_cor_norm_idx << endl<< endl; 
                    else
                        levWrite(std::cout,mg_level,str) << '\t' << "residual eval only - no correction" << endl<< endl; 
                }
            }
        } // nu2 iters
   } // if not a V-bottom            
}

void
HeatTransfer::mcdd_update(Real time,
                          Real dt)
{
    BL_ASSERT(do_mcdd);

    DDOp::DD_ApForTorRH whichApp = (mcdd_advance_temp==1 ? DDOp::DD_Temp : DDOp::DD_RhoH);
    std::string vgrp("mcdd"); // For verbose debugging info

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mcdd update for "
                  << (whichApp==DDOp::DD_Temp ? "Temp" : "RhoH") 
                  << " and RhoY" << '\n';

    showMF(vgrp,*aofs,"aofs_mcdd_update");

    // Get advection updates into new-time state
    scalar_advection_update(dt, Density, Density);
    scalar_advection_update(dt, first_spec, last_spec);
    scalar_advection_update(dt, RhoH, RhoH);
    scalar_advection_update(dt, Temp, Temp);

    if (hack_nospecdiff==1)
        return;

    //
    //   If whichApp says do T, here's the plan.  The temperature
    //     equation is discretized as follows (HS(t)=heat source at t,
    //      such as Div(Qrad)):
    //
    //   (rho.cp)^{nph} ( Tnew - Told + dt.u.grad(T) ) =
    //         theta .dt.( (1/V) Div(QT) - Fi.cpi.Grad(T) )^{np1}
    //      (1-theta).dt.( (1/V) Div(QT) - Fi.cpi.Grad(T) )^{n}
    //         theta.dt.HS(tnew) + (1-theta).dt.HS(told)
    //
    //   Here, QT is the heat flux without the Fi.hi contributions.
    //   Moving things around, we get
    //
    //                  Tnew - theta.dt.LT(tnew) = Rhs
    //
    //   where Rhs = Told - dt.u.Grad(T) 
    //        + [ (1-theta)( LT(told) + HS(told) ) + theta.HS(tnew) ].dt/(rho.cp)^{nph}
    //   
    //   Here, LT(t) is what DDOp.apply returns if whichApp is for T.  Also,
    //     S_new = S_old - dt*aofs, where aofs(T) usully is u.grad(T)
    //   Really, who knows what advection update might be adding in?
    //     Whatever it did, we'll treat it as the "advection update"
    //
    //   If whichApp says to do RhoH instead, there is very little change.
    //
    //      ( RhoHnew - RhoHold + dt.Div(u.RhoH) ) =
    //         theta .dt.( (1/V) Div(Q) )^{np1}
    //      (1-theta).dt.( (1/V) Div(Q) )^{n}
    //         theta.dt.HS(tnew) + (1-theta).dt.HS(told)
    //   

    const int nGrow = 0;
    MultiFab Rhs(grids,nspecies+1,nGrow); //stack H on Y's
    const int dCompY = 0;
    const int dCompH = nspecies;

    // Build (Y,T,H,rho) first with old data to evaluate Rhs
    MultiFab& S_old = get_old_data(State_Type);
    showMF(vgrp,S_old,"S_old_mcdd_update");

    const int nGrowOp = 1;
    MultiFab YTHR(grids,nspecies+3,nGrowOp);
    MultiFab T1;
    if (whichApp == DDOp::DD_Temp) {
        T1.define(grids,1,0,Fab_allocate); // Holds cpold, then rhoCpInv^{nph}
    }
    int dCompSY = 0;
    int dCompST = dCompSY + nspecies;
    int dCompSH = dCompST + 1;
    int dCompSR = dCompSH + 1;
    for (MFIter mfi(YTHR); mfi.isValid(); ++mfi)
    {
        FArrayBox&    ythr = YTHR[mfi];
        const FArrayBox& s = S_old[mfi];
        const Box& box = mfi.validbox();

        ythr.copy(s,box,Density,box,dCompSR,1);        
        ythr.copy(s,box,Temp,box,dCompST,1);
        ythr.copy(s,box,first_spec,box,dCompSY,nspecies);
        for (int i=0; i<nspecies; ++i)
            ythr.divide(ythr,box,dCompSR,dCompSY+i,1);
        ythr.copy(s,box,RhoH,box,dCompSH,1);
        ythr.divide(ythr,box,dCompSR,dCompSH,1);

        // While we have a (Y,T)_old, get cp_old
        if (whichApp == DDOp::DD_Temp)
        {
            const int sCompCp = 0;
            getChemSolve().getCpmixGivenTY(T1[mfi],ythr,ythr,box,dCompST,dCompSY,sCompCp);
        }
    }
    showMF(vgrp,YTHR,"YTHR_old_mcdd_update");
    if (whichApp == DDOp::DD_Temp)
        showMF(vgrp,T1,"cp_old_mcdd_update");

    PArray<MultiFab> fluxn(BL_SPACEDIM,PArrayManage);
    for (int dir=0; dir<BL_SPACEDIM; ++dir)
        fluxn.set(dir,new MultiFab(BoxArray(grids).surroundingNodes(dir),
                                   nspecies+1,0));

    // Load boundary data in operator for old time, compute L(told)
    fill_mcdd_boundary_data(time);
    MCDDOp.setCoefficients(YTHR,dCompST,YTHR,dCompSY);
    showMCDD("mcddOp",MCDDOp,"MCDD_old");
    int level = 0;
    bool getAlpha = false;
    MultiFab* alpha = 0;
    MCDDOp.applyOp(Rhs,YTHR,fluxn,whichApp,level,getAlpha,alpha);
    showMF(vgrp,Rhs,"L_of_YTHR_old_mcdd_update");

    add_heat_sources(Rhs,dCompH,time,0,1);
    Rhs.mult(1-be_cn_theta);

    add_heat_sources(Rhs,dCompH,time+dt,0,be_cn_theta); // Note: based on adv-predicted state
    showMF(vgrp,Rhs,"Rhs_old_stuff_mcdd_update");

    // Build YTHR_new, and complete the build of Rhs
    MultiFab& S_new = get_new_data(State_Type);
    showMF(vgrp,S_new,"S_new_mcdd_update");

    FArrayBox cp_new, rho_half;
    for (MFIter mfi(YTHR); mfi.isValid(); ++mfi)
    {
        const Box&     box = mfi.validbox();
        FArrayBox&    ythr = YTHR[mfi];
        const FArrayBox& s = S_new[mfi];

        ythr.copy(s,box,Density,box,dCompSR,1);        
        ythr.copy(s,box,Temp,box,dCompST,1);
        ythr.copy(s,box,first_spec,box,dCompSY,nspecies);
        for (int i=0; i<nspecies; ++i)
            ythr.divide(ythr,box,dCompSR,dCompSY+i,1);
        ythr.copy(s,box,RhoH,box,dCompSH,1);
        ythr.divide(ythr,box,dCompSR,dCompSH,1);

        if (whichApp == DDOp::DD_Temp)
        {        
            cp_new.resize(box,1);
            const int sCompCp = 0;
            getChemSolve().getCpmixGivenTY(cp_new,ythr,ythr,box,dCompST,dCompSY,sCompCp);
            T1[mfi].plus(cp_new,0,0,1);  // T1 currently cp_old + cp_new
            
            rho_half.resize(box,1);
            rho_half.copy(S_old[mfi],Density,0,1);
            rho_half.plus(S_new[mfi],Density,0,1);
            
            T1[mfi].mult(rho_half,0,0,1);
            T1[mfi].invert(4); // Note: 2 from cp average, and 2 from rho average

            Rhs[mfi].mult(T1[mfi],0,dCompH,1);
        }

        Rhs[mfi].mult(dt);

        // Add in the S_new = S_old - dt.aofs term
        Rhs[mfi].plus(S_new[mfi],first_spec,dCompY,nspecies);
        if (whichApp==DDOp::DD_Temp) {
            Rhs[mfi].plus(S_new[mfi],Temp,dCompH,1);
        } else {
            Rhs[mfi].plus(S_new[mfi],RhoH,dCompH,1);
        }
    }

    showMF(vgrp,YTHR,"YTHR_new_mcdd_update");
    showMF(vgrp,Rhs,"Rhs_mcdd_update");
    if (whichApp == DDOp::DD_Temp)
        showMF(vgrp,T1,"rhoCpInv_mcdd_update");

    PArray<MultiFab> fluxnp1(BL_SPACEDIM,PArrayManage);
    for (int dir=0; dir<BL_SPACEDIM; ++dir)
        fluxnp1.set(dir,new MultiFab(BoxArray(grids).surroundingNodes(dir),
                                     nspecies+1,nGrow)); //stack H on Y's

    // Fill the boundary data for the solve (all the levels)
    fill_mcdd_boundary_data(time+dt);
    MCDDOp.setCoefficients(YTHR,dCompST,YTHR,dCompSY);
    showMCDD("mcddOp",MCDDOp,"MCDD_new");

    MCDD_MGParams p(nspecies+1);
    p.stalled_tol = mcdd_stalled_tol;
    p.res_nu1_rtol = mcdd_res_nu1_rtol;
    p.res_nu2_rtol = mcdd_res_nu2_rtol;
    p.res_redux_tol = mcdd_res_redux_tol;
    p.res_abs_tol = mcdd_res_abs_tol;
    p.mg_level = 0;
    p.gamma = 1; // FIXME: Want as input?
    p.num_coarser = DDOp::numLevels() - 1;
    if (p.num_coarser<1  &&  mcdd_nub>0) {
        p.nu1 = mcdd_nub;
        p.nu2 = 0;
    } else {
        p.nu1 = mcdd_nu1[p.mg_level];
        p.nu2 = mcdd_nu2[p.mg_level];
    }
    p.status = HT_InProgress;

    bool update_mcdd_coeffs = true;
    //bool update_mcdd_coeffs = false;
    for (p.iter=0; (p.iter<mcdd_numcycles) && (p.status==HT_InProgress); ++p.iter)
    {
        mcdd_fas_cycle(p,YTHR,Rhs,T1,fluxnp1,time,dt,whichApp);

        if (ParallelDescriptor::IOProcessor() && (mcdd_verbose>0) ) {
            std::string str("mcdd - (iter,res,cor): (");
            levWrite(std::cout,-1,str)
                << p.iter << ", "
                << p.maxRes_norm << ", "
                << p.maxCor_norm << ") " << std::endl;
        }

        if (update_mcdd_coeffs && p.status==HT_InProgress)
        {
            sync_H_T(YTHR,getChemSolve(),whichApp,dCompSY,dCompST,dCompSH,mcdd_verbose,0," (outer loop) ");

            if (ParallelDescriptor::IOProcessor()) {
                std::string str("mcdd - updating transport coeffs ...");
                levWrite(std::cout,-1,str) << endl;
            }

            MCDDOp.setCoefficients(YTHR,dCompST,YTHR,dCompSY);
            //showMCDD("mcddOp",MCDDOp,"MCDD_new");
        }
    }

    if (ParallelDescriptor::IOProcessor()) {
        if (p.status == HT_InProgress) {
            std::cout << "**************** mcdd solver failed after " << mcdd_numcycles << " V-cycles!!" << std::endl;

            cout << "Residual redux (R,Ri,R/Ri): " << endl;
            for (int i=0; i<nspecies+1; ++i) {
                cout << i << " "
                     << p.maxRes[i] << " "
                     << p.maxRes_initial[i] << " "
                     << p.maxRes[i]/p.maxRes_initial[i] << endl;
            }
        }
    }

    // Put solution back into state (after scaling by rhonew)
    for (MFIter mfi(YTHR); mfi.isValid(); ++mfi)
    {
        FArrayBox& s = S_new[mfi];
        const FArrayBox& ythr = YTHR[mfi];
        const Box& box = mfi.validbox();

        s.copy(ythr,box,dCompSR,box,Density,1);        
        s.copy(ythr,box,dCompST,box,Temp,1);
        s.copy(ythr,box,dCompSY,box,first_spec,nspecies);
        for (int i=0; i<nspecies; ++i)
            s.mult(s,box,Density,first_spec+i,1);
        s.copy(ythr,box,dCompSH,box,RhoH,1);
        s.mult(s,box,Density,RhoH,1);
    }

    if (do_reflux)
    {
        //
        // TODO - integrate in new CrseInit() call.
        //

        // Build the total flux (pieces from n and np1) for RY and RH
        FArrayBox fluxtot, fluxtmp;
        for (int d = 0; d < BL_SPACEDIM; d++)
        {
            for (MFIter mfi(fluxn[d]); mfi.isValid(); ++mfi)
            {
                const FArrayBox& fn = fluxn[d][mfi];
                const FArrayBox& fnp1 = fluxnp1[d][mfi];
                const Box& ebox = fn.box();

                fluxtot.resize(ebox,nspecies+1);
                fluxtmp.resize(ebox,nspecies+1);
                fluxtot.copy(fn,ebox,0,ebox,0,nspecies+1);
                fluxtot.mult(1.0-be_cn_theta);
                fluxtmp.copy(fnp1,ebox,0,ebox,0,nspecies+1);
                fluxtmp.mult(be_cn_theta);
                fluxtot.plus(fluxtmp,ebox,0,0,nspecies+1);
                if (level < parent->finestLevel())
                {
                    // Note: The following inits do not trash eachother
                    //  since the component ranges don't overlap...
                    getLevel(level+1).getViscFluxReg().CrseInit(fluxtot,ebox,
                                                                d,dCompY,first_spec,
                                                                nspecies,-dt);
                    getLevel(level+1).getViscFluxReg().CrseInit(fluxtot,ebox,
                                                                d,dCompH,RhoH,
                                                                1,-dt);
                }
                if (level > 0)
                {
                    getViscFluxReg().FineAdd(fluxtot,d,mfi.index(),
                                             dCompY,first_spec,nspecies,dt);
                    getViscFluxReg().FineAdd(fluxtot,d,mfi.index(),
                                             dCompH,RhoH,1,dt);
                }
            }
        }
        if (level < parent->finestLevel())
            getLevel(level+1).getViscFluxReg().CrseInitFinish();
    }
}

static
int HY_to_Temp_DoIt(FArrayBox&       Tfab,
                    const FArrayBox& Hfab,
                    const FArrayBox& Yfab,
                    const Box&       box,
                    int              sCompH,
                    int              sCompY,
                    int              dCompT,
                    Real             htt_hmixTYP,
                    ChemDriver&      cd)
{
    if (htt_hmixTYP <= 0)
    {
        htt_hmixTYP = Hfab.norm(0,sCompH,1);
        if (ParallelDescriptor::IOProcessor())
            std::cout << "setting htt_hmixTYP = " << htt_hmixTYP << '\n';        
    }

    const Real eps = cd.getHtoTerrMAX();
    Real errMAX = eps*htt_hmixTYP;

    int iters = cd.getTGivenHY(Tfab,Yfab,Hfab,box,sCompH,sCompY,dCompT,errMAX);

    if (iters < 0)
        BoxLib::Error("HeatTransfer::RhoH_to_Temp(fab): error in H->T");

    return iters;
}

void
HeatTransfer::compute_differential_diffusion_fluxes (const Real& time,
                                                     const Real& dt)
{
    const TimeLevel whichTime = which_time(State_Type,time);
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    
    MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;

    int rhoD_rhs_comp = 0; // Starting slot for beta, rhs for diffusion calls

    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! skipping spec diffusion " << '\n';
        
        for (int d = 0; d < BL_SPACEDIM; ++d)
        {
            (*flux[d]).setVal(0,rhoD_rhs_comp,nspecies+1);
        }

        return;
    }

    MultiFab* rho_half = 0; // Never need alpha for RhoY
    MultiFab* alpha = 0;
    int alphaComp = 0;
    const Real a = 0;
    const Real b = 1;
    const int rho_flag = 2;
    MultiFab **beta = 0;
    diffusion->allocFluxBoxesLevel(beta,0,nspecies+2); // Species and heat fluxes
    getDiffusivity(beta, time, first_spec, 0, nspecies+1);
    getDiffusivity(beta, time, Temp, nspecies+1, 1);
    MultiFab& S = get_data(State_Type,time);

    // Fill grow cells in the state for (Rho,RhoY,RhoH,T)
    // For Dirichlet physical boundary grow cells, this data will live on the cell face, otherwise
    // it will live at cell centers
    int nGrow = 1;
    BL_ASSERT(S.nGrow()>=nGrow);
    FillPatchIterator Yfpi(*this,S,nGrow,time,State_Type,Density,nspecies+2);
    FillPatchIterator Tfpi(*this,S,nGrow,time,State_Type,Temp,1);
    for ( ; Yfpi.isValid() && Tfpi.isValid(); ++Yfpi, ++Tfpi)
    {
        const Box& vbox = Yfpi.validbox();
        FArrayBox& fab = S[Tfpi];
        BoxList gcells = BoxLib::boxDiff(Box(vbox).grow(nGrow),vbox);
        for (BoxList::const_iterator it = gcells.begin(), end = gcells.end(); it != end; ++it)
        {
            const Box& gbox = *it;
            fab.copy(Yfpi(),gbox,0,gbox,Density,nspecies+2);
            fab.copy(Tfpi(),gbox,0,gbox,Temp,1);
        }
    }
    showMF("dd",S,"dd_rsT_fp",level);
    //
    // Create and fill (Rho,RhoY,RhoH,T) at coarser level
    //    
    MultiFab* S_crse = 0;
    MultiFab* Phi_crse = 0; 
    if (level > 0) 
    {
      HeatTransfer* coarser = (HeatTransfer*) &(parent->getLevel(level-1));


      BoxArray crse_grids(coarser->grids);
      Phi_crse = new MultiFab(crse_grids, 1, 2);
      S_crse = new MultiFab(crse_grids,S.nComp(),0,Fab_allocate);

      for (FillPatchIterator S_fpi(*coarser,*S_crse,0,time,State_Type,0,S.nComp());
	   S_fpi.isValid(); ++S_fpi)
	{
	  (*S_crse)[S_fpi.index()].copy(S_fpi(),0,0,S.nComp());
	}
    }

    MultiFab Phi(grids,1,nGrow);
    int phiComp = 0;

    for (int sigma = 0; sigma < nspecies+1; ++sigma)
    {
        const int state_ind = first_spec + sigma;
        ViscBndry visc_bndry;

        for (MFIter mfi(S); mfi.isValid(); ++mfi) {
            const Box& box = Phi[mfi].box();
            Phi[mfi].copy(S[mfi],box,state_ind,box,0,1);
            Phi[mfi].divide(S[mfi],box,Density,0,1);
        }

        if (level > 0) 
        {
            MultiFab::Copy(*Phi_crse,*S_crse,state_ind,0,1,0);
            for (MFIter mfi(*S_crse); mfi.isValid(); ++mfi)
	      {
		(*Phi_crse)[mfi].divide((*S_crse)[mfi],(*Phi_crse)[mfi].box(),Density,0,1);
	      }
	    Phi_crse->FillBoundary();
	    getLevel(level-1).geom.FillPeriodicBoundary(*Phi_crse);
        }

	diffusion->getBndryDataGivenS(visc_bndry,Phi,*Phi_crse,state_ind,0,1);

        bool bndry_already_filled = true;            
        int betaComp = sigma;
        Real* rhsscale = 0;
        ABecLaplacian* visc_op = diffusion->getViscOp(state_ind,a,b,time,visc_bndry,
                                                      rho_half,rho_flag,rhsscale,
                                                      beta,betaComp,alpha,alphaComp,bndry_already_filled);
        visc_op->maxOrder(diffusion->maxOrder());
        bool do_applyBC = true;
        visc_op->compFlux(D_DECL(*flux[0],*flux[1],*flux[2]),Phi,LinOp::Inhomogeneous_BC,do_applyBC,phiComp,sigma);
        delete visc_op;
    }

    if (Phi_crse)
      delete Phi_crse;

    if (S_crse)
      delete S_crse;

    // Remove scaling left in fluxes from solve
    for (int d=0; d < BL_SPACEDIM; ++d) {
        flux[d]->mult(b/geom.CellSize()[d],0,nspecies+1);
    }
    //
    // Modify update/fluxes to preserve flux sum = 0, build heat flux and temperature source terms
    //
    bool grow_cells_already_filled = true;

    // conservatively correct Gamma_m
    adjust_spec_diffusion_fluxes(time,grow_cells_already_filled);

    //
    // AJN FLUXREG
    // We have just explicitly computed "D" for Y_m and h given an input state.
    // Update VISCOUS flux registers as follows.
    // If we are calling this in the predictor:
    //   COPY (1/2)*Gamma_m^n and (1/2)*lambda^n/cp^n grad h^n to flux registers
    // If we are calling this in the corrector AND updateFluxReg=T:
    //   SUBTRACT (1/2)*Gamma_m^(k) and (1/2)*lambda^(k)/cp^(k) grad h^(k)
    if ( (do_reflux && is_predictor) ||
         (do_reflux && !is_predictor && updateFluxReg) )
    {
      const Real theta = (is_predictor) ? 0.5 : -0.5;

      for (int d = 0; d < BL_SPACEDIM; d++)
      {
	for (MFIter fmfi(*SpecDiffusionFluxn[d]); fmfi.isValid(); ++fmfi)
	{
	  if (level > 0)
          {
	    if (is_predictor)
	      getViscFluxReg().FineAdd((*SpecDiffusionFluxn[d])[fmfi],d,fmfi.index(),0,first_spec,nspecies+1,theta*dt);
	    else
	      getViscFluxReg().FineAdd((*SpecDiffusionFluxnp1[d])[fmfi],d,fmfi.index(),0,first_spec,nspecies+1,theta*dt);
	  }
	}

	if (level < parent->finestLevel())
        {
	  if (is_predictor)
	    getLevel(level+1).getViscFluxReg().CrseInit((*SpecDiffusionFluxn[d]),d,0,first_spec,nspecies+1,-theta*dt,FluxRegister::COPY);
	  else
	    getLevel(level+1).getViscFluxReg().CrseInit((*SpecDiffusionFluxnp1[d]),d,0,first_spec,nspecies+1,-theta*dt,FluxRegister::ADD);
	}
      }
    }

    // compute lambda grad T (for temperature and divu)
    // compute sum_m (Gamma_m + lambda/cp grad Y) (for enthalpy)
    // compute sum_m Gamma_m dot grad h_m (for divu)
    compute_enthalpy_fluxes(time,beta,grow_cells_already_filled);

    //
    // AJN FLUXREG
    // We have just explicitly computed "DD" given an input state.
    // Update ADVECTIVE flux registers as follows.
    // If we are calling this in the predictor:
    //   ADD (1/2)*h_m^n(Gamma_m^n-lambda^n/cp^n grad Y_m^n)
    // If we are calling this in the corrector AND updateFluxReg=T:
    //   ADD (1/2)*h_m^(k)(Gamma_m^(k)-lambda^(k)/cp^(k) grad Y_m^(k))
    if ( (do_reflux && is_predictor) ||
         (do_reflux && !is_predictor && updateFluxReg) )
    {
      for (int d = 0; d < BL_SPACEDIM; d++)
      {
	for (MFIter fmfi(*SpecDiffusionFluxn[d]); fmfi.isValid(); ++fmfi)
	{
	  if (level > 0)
          {
	    if (is_predictor)
	      advflux_reg->FineAdd((*SpecDiffusionFluxn[d])[fmfi],  d,fmfi.index(),nspecies+1,RhoH,1,0.5*dt);
	    else
	      advflux_reg->FineAdd((*SpecDiffusionFluxnp1[d])[fmfi],d,fmfi.index(),nspecies+1,RhoH,1,0.5*dt);
	  }
	}

	if (level < parent->finestLevel())
        {
	  if (is_predictor)
	    getAdvFluxReg(level+1).CrseInit((*SpecDiffusionFluxn[d])  ,d,nspecies+1,RhoH,1,-0.5*dt,FluxRegister::ADD);
	  else
	    getAdvFluxReg(level+1).CrseInit((*SpecDiffusionFluxnp1[d]),d,nspecies+1,RhoH,1,-0.5*dt,FluxRegister::ADD);
	}
      }
    }

    diffusion->removeFluxBoxesLevel(beta);
}

void
HeatTransfer::scalar_advection_update (Real dt,
                                       int  first_scalar,
                                       int  last_scalar)
{
    // Careful: If here, the sign of aofs is flipped (wrt the usual NS treatment)
  FArrayBox tmp;
  MultiFab& S_new = get_new_data(State_Type);
  const MultiFab& S_old = get_old_data(State_Type);
  
  for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
      const Box& box   = mfi.validbox();
      const int nc = last_scalar - first_scalar + 1; 
      FArrayBox& snew = S_new[mfi];
      
      snew.copy( (*aofs)[mfi],box,first_scalar,box,first_scalar,nc);
      snew.mult(dt,first_scalar,nc);
      snew.plus(S_old[mfi],first_scalar,first_scalar,nc);            
    }
}

void
HeatTransfer::flux_divergence (MultiFab&        fdiv,
                               int              fdivComp,
                               const MultiFab* const* f,
                               int              fluxComp,
                               int              nComp,
                               Real             scale) const
{
    BL_ASSERT(fdiv.nComp() >= fdivComp+nComp);

    FArrayBox volume;
    for (MFIter mfi(fdiv); mfi.isValid(); ++mfi)
    {
	int        i = mfi.index();
	const Box& box = mfi.validbox();
        FArrayBox& fab = fdiv[mfi];
        geom.GetVolume(volume,grids,i,GEOM_GROW);

	FORT_FLUXDIV(box.loVect(), box.hiVect(),
                     fab.dataPtr(fdivComp), ARLIM(fab.loVect()), ARLIM(fab.hiVect()),
                     (*f[0])[i].dataPtr(fluxComp), ARLIM((*f[0])[i].loVect()),   ARLIM((*f[0])[i].hiVect()),
                     (*f[1])[i].dataPtr(fluxComp), ARLIM((*f[1])[i].loVect()),   ARLIM((*f[1])[i].hiVect()),
#if BL_SPACEDIM == 3
                     (*f[2])[i].dataPtr(fluxComp), ARLIM((*f[2])[i].loVect()),   ARLIM((*f[2])[i].hiVect()),
#endif
                     volume.dataPtr(),             ARLIM(volume.loVect()),       ARLIM(volume.hiVect()),
                     &nComp, &scale);
    }	
}

void
HeatTransfer::compute_differential_diffusion_terms (MultiFab& D,
                                                    MultiFab& DD,
                                                    Real      time,
                                                    Real      dt)
{
    // 
    // Sets vt for species, RhoH and Temp together
    // Uses state at time to explicitly compute fluxes, and resets internal
    //  data for fluxes, etc
    //
    int nComp = nspecies + 2;
    BL_ASSERT(D.boxArray() == grids);
    BL_ASSERT(D.nComp() >= nComp);// spec+heat+conduction

    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! zeroing diffusion terms " << '\n';
        D.setVal(0.0,0,nComp);
        return;            
    }

    const TimeLevel whichTime = which_time(State_Type,time);
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);    
    MultiFab& sumSpecFluxDotGradH = (whichTime == AmrOldTime) ? sumSpecFluxDotGradHn : sumSpecFluxDotGradHnp1;
    MultiFab* const * flux = (whichTime == AmrOldTime) ? SpecDiffusionFluxn : SpecDiffusionFluxnp1;
    //
    // Compute/adjust species fluxes/heat flux/conduction, save in class data
    // Then compute -Div(Fi), -Div((lambda/cp).Grad(h) + Fi.(Lei-1)) + heating,
    //  -Div(lambda.Grad(T)) + heating + sum(Fi.Grad(Hi)) 
    //
    compute_differential_diffusion_fluxes(time,dt);

    // compute div Gamma_m for species AND div lambda/cp grad h for enthalpy
    flux_divergence(D,0,flux,0,nspecies+1,-1);

    // compute div sum_m h_m (Gamma_m + lambda/cp grad Y_m), a.k.a. the "diffdiff" terms
    flux_divergence(DD,0,flux,nspecies+1,1,-1);

    // compute div lambda grad T for temperature
    flux_divergence(D,nspecies+1,flux,nspecies+2,1,-1);

    // add sum_m Gamma_m dot grad h_m to D for temperature
    MultiFab::Add(D,sumSpecFluxDotGradH,0,nspecies+1,1,0);

    if (D.nGrow() > 0 && DD.nGrow() > 0)
      {
	BL_ASSERT(D.nGrow() == DD.nGrow());

	const int nc = nspecies+2;
	const int ncDD = 1;
	for (MFIter mfi(D); mfi.isValid(); ++mfi)
	  {
	    FArrayBox& Dfab = D[mfi];
	    FArrayBox& DDfab = DD[mfi];
	    const Box& box = mfi.validbox();
	    FORT_VISCEXTRAP(Dfab.dataPtr(),ARLIM(Dfab.loVect()),ARLIM(Dfab.hiVect()),
			    box.loVect(),box.hiVect(),&nc);
	    FORT_VISCEXTRAP(DDfab.dataPtr(),ARLIM(DDfab.loVect()),ARLIM(DDfab.hiVect()),
			    box.loVect(),box.hiVect(),&ncDD);
	  }
	D.FillBoundary(0,nc);
	geom.FillPeriodicBoundary(D,0,nc,true);
	DD.FillBoundary(0,1);
	geom.FillPeriodicBoundary(DD,0,1,true);
      }
}

void
HeatTransfer::add_heat_sources(MultiFab& sum,
                               int       dComp,
                               Real      time,
                               int       nGrow,
                               Real      scale)
{
    //
    // - div q rad
    //
    if (do_OT_radiation || do_heat_sink)
    {
        BL_ASSERT(sum.boxArray() == grids);
        MultiFab dqrad(grids,1,nGrow);
        compute_OT_radloss(time,nGrow,dqrad);
        for (MFIter mfi(sum); mfi.isValid(); ++mfi)
        {
            FArrayBox& sumFab = sum[mfi];
            FArrayBox& Qfab = dqrad[mfi];
            if (scale != 1)
                Qfab.mult(scale);
            sumFab.minus(Qfab,mfi.validbox(),0,dComp,1);
        }
    }    
}

void
HeatTransfer::set_rho_to_species_sum (MultiFab& S,
                                      int       strtcomp, 
                                      int       nghost_in,
                                      int       minzero)
{
    set_rho_to_species_sum(S, strtcomp, S, strtcomp, nghost_in, minzero);
}

//
// This function
//       sets the density component in S_out to the sum of 
//       the species components in S_in
// if minzero = 1, the species components are first "max-ed" w/ zero
//  thus changing S_in values    
//
// s_in_strt is the state component corresponding to the
// 0-th component of S_in. It is otherwise assumed that 
// that the components of S_in "align" with those in S_out.
//
void
HeatTransfer::set_rho_to_species_sum (MultiFab& S_in, 
                                      int       s_in_strt,
                                      MultiFab& S_out,
                                      int       s_out_strt,  
                                      int       nghost_in,
                                      int       minzero)

{
    const BoxArray& sgrids = S_in.boxArray();

    BL_ASSERT(sgrids == S_out.boxArray());

    const int s_first_spec = first_spec - s_in_strt;
    const int s_last_spec  = last_spec  - s_in_strt;
    const int s_num_spec   = last_spec - first_spec + 1;
    const int s_density    = Density-s_out_strt;
    const int nghost       = std::min(S_in.nGrow(), std::min(S_out.nGrow(), nghost_in));

    if (minzero)
    {
        for (MFIter mfi(S_in); mfi.isValid(); ++mfi)
        {
            const int i   = mfi.index();
            Box       box = BoxLib::grow(sgrids[i], nghost);
            FabMinMax(S_in[i], box, 0.0, Real_MAX, s_first_spec, s_num_spec);
        }
    }

    BL_ASSERT(s_density >= 0);

    S_out.setVal(0, s_density, 1, nghost);

    for (MFIter mfi(S_in); mfi.isValid(); ++mfi)
    {
        const int i   = mfi.index();
        Box       box = BoxLib::grow(sgrids[i],nghost);

        for (int spec = s_first_spec; spec <= s_last_spec; spec++)
        {
            S_out[i].plus(S_in[i],box,spec,s_density,1);
        }
    }

    make_rho_curr_time();
}

void
HeatTransfer::temperature_stats (MultiFab& S)
{
    if (verbose)
    {
        //
        // Calculate some minimums and maximums.
        //
        Real tdhmin[3] = { 1.0e30, 1.0e30, 1.0e30};
        Real tdhmax[3] = {-1.0e30,-1.0e30,-1.0e30};

        for (MFIter S_mfi(S); S_mfi.isValid(); ++S_mfi)
        {
            const Box& bx = S_mfi.validbox();

            tdhmin[0] = std::min(tdhmin[0],S[S_mfi].min(bx,Temp));
            tdhmin[1] = std::min(tdhmin[1],S[S_mfi].min(bx,Density));
            tdhmin[2] = std::min(tdhmin[2],S[S_mfi].min(bx,RhoH));

            tdhmax[0] = std::max(tdhmax[0],S[S_mfi].max(bx,Temp));
            tdhmax[1] = std::max(tdhmax[1],S[S_mfi].max(bx,Density));
            tdhmax[2] = std::max(tdhmax[2],S[S_mfi].max(bx,RhoH));
        }

        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceRealMin(tdhmin,3,IOProc);
        ParallelDescriptor::ReduceRealMax(tdhmax,3,IOProc);

        FArrayBox   Y, tmp;
        bool        aNegY = false;
        Array<Real> minY(nspecies,1.e30);

        for (MFIter S_mfi(S); S_mfi.isValid(); ++S_mfi)
        {
            Y.resize(S_mfi.validbox(),1);

            tmp.resize(S_mfi.validbox(),1);
            tmp.copy(S[S_mfi],Density,0,1);
            tmp.invert(1);

            for (int i = 0; i < nspecies; ++i)
            {
                Y.copy(S[S_mfi],first_spec+i,0,1);
                Y.mult(tmp,0,0,1);
                minY[i] = std::min(minY[i],Y.min(0));
            }
        }

        ParallelDescriptor::ReduceRealMin(minY.dataPtr(),nspecies,IOProc);

        for (int i = 0; i < nspecies; ++i)
        {
            if (minY[i] < 0) aNegY = true;
        }

        if (ParallelDescriptor::IOProcessor())
        {
            std::cout << "  Min,max temp = " << tdhmin[0] << ", " << tdhmax[0] << '\n';
            std::cout << "  Min,max rho  = " << tdhmin[1] << ", " << tdhmax[1] << '\n';
            std::cout << "  Min,max rhoh = " << tdhmin[2] << ", " << tdhmax[2] << '\n';
            if (aNegY)
            {
                const Array<std::string>& names = HeatTransfer::getChemSolve().speciesNames();
                std::cout << "  Species w/min < 0: ";
                for (int i = 0; i < nspecies; ++i)
                    if (minY[i] < 0)
                        std::cout << "Y(" << names[i] << ") [" << minY[i] << "]  ";
                std::cout << '\n';
            }
        }
    }
}

void
HeatTransfer::compute_rhoRT (const MultiFab& S,
                             MultiFab&       p, 
                             int             pComp,
                             const MultiFab* temp)
{
    BL_ASSERT(pComp<p.nComp());

    const BoxArray& sgrids = S.boxArray();
    int nCompY = last_spec - first_spec + 1;

    FArrayBox state, tmp;

    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const int  i      = mfi.index();
	const Box& box    = sgrids[i];
	const int  sCompR = 0;
	const int  sCompT = 1;
	const int  sCompY = 2;
	
        state.resize(box,nCompY+2);
	BL_ASSERT(S[mfi].box().contains(box));
	state.copy(S[mfi],box,Density,box,sCompR,1);

	if (temp)
	{
	    BL_ASSERT(temp->boxArray()[i].contains(box));
	    state.copy((*temp)[i],box,0,box,sCompT,1);
	}
        else
        {
	    state.copy(S[mfi],box,Temp,box,sCompT,1);
	}
	state.copy(S[mfi],box,first_spec,box,sCompY,nCompY);

        tmp.resize(box,1);
        tmp.copy(state,sCompR,0,1);
        tmp.invert(1);

        for (int k = 0; k < nCompY; k++)
            state.mult(tmp,0,sCompY+k,1);

	getChemSolve().getPGivenRTY(p[i],state,state,state,box,sCompR,sCompT,sCompY,pComp);
    }
}
			   
//
// Setup for the advance function.
//

#ifndef NDEBUG
#if defined(BL_OSF1)
#if defined(BL_USE_DOUBLE)
const Real BL_BOGUS      = DBL_QNAN;
#else
const Real BL_BOGUS      = FLT_QNAN;
#endif
#else
const Real BL_BOGUS      = 1.e30;
#endif
#endif

void
HeatTransfer::set_htt_hmixTYP ()
{
    const int finest_level = parent->finestLevel();

    // set typical value for hmix, needed for TfromHY solves if not provided explicitly
    if (typical_values[RhoH] < 0) {
        htt_hmixTYP = -1;
        for (int k = 0; k <= finest_level; k++)
        {
            AmrLevel& ht = getLevel(k);
            const MultiFab& S = ht.get_new_data(State_Type);
            const BoxArray& ba = ht.boxArray();
            MultiFab hmix(ba,1,0,Fab_allocate);
            MultiFab::Copy(hmix,S,RhoH,0,1,0);
            MultiFab::Divide(hmix,S,Density,0,1,0);
            if (k!=finest_level) 
            {
                AmrLevel& htf = getLevel(k+1);
                BoxArray baf = BoxArray(htf.boxArray()).coarsen(parent->refRatio(k));
                for (MFIter mfi(hmix); mfi.isValid(); ++mfi)
                {
                    std::vector< std::pair<int,Box> > isects = baf.intersections(ba[mfi.index()]);
                    for (int i = 0; i < isects.size(); i++)
                    {
                        hmix[mfi].setVal(0,isects[i].second,0,1);
                    }
                }
            }
            htt_hmixTYP = std::max(htt_hmixTYP,hmix.norm0(0));
        }
        ParallelDescriptor::ReduceRealMax(htt_hmixTYP);
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "setting htt_hmixTYP(via domain scan) = " << htt_hmixTYP << '\n';
    } else {
        htt_hmixTYP = typical_values[RhoH];
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "setting htt_hmixTYP(from user input) = " << htt_hmixTYP << '\n';
    }
}

void
HeatTransfer::advance_setup (Real time,
                             Real dt,
                             int  iteration,
                             int  ncycle)
{
    NavierStokes::advance_setup(time, dt, iteration, ncycle);
    //
    // Make sure the new state has values so that c-n works
    // in the predictor (of the predictor)--rbp.
    //
    for (int k = 0; k < num_state_type; k++)
    {
        MultiFab& nstate = get_new_data(k);
        MultiFab& ostate = get_old_data(k);

        MultiFab::Copy(nstate,ostate,0,0,nstate.nComp(),0);
    }
    if (level == 0)
        set_htt_hmixTYP();

    make_rho_curr_time();

#ifndef NDEBUG
    aux_boundary_data_old.setVal(BL_BOGUS);
    aux_boundary_data_new.setVal(BL_BOGUS);
#endif
}

void
HeatTransfer::setThermoPress(Real time)
{
    const TimeLevel whichTime = which_time(State_Type,time);
    
    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);
    
    MultiFab& S = (whichTime == AmrOldTime) ? get_old_data(State_Type) : get_new_data(State_Type);
    
    const int pComp = (have_rhort ? RhoRT : Trac);

    compute_rhoRT (S,S,pComp);
}

#if 1
static int hyp_grow = 3; // ick!  
Real
HeatTransfer::predict_velocity (Real  dt,
                                Real& comp_cfl)
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... predict edge velocities\n";
    //
    // Get simulation parameters.
    //
    const int   nComp          = BL_SPACEDIM;
    const Real* dx             = geom.CellSize();
    const Real  prev_time      = state[State_Type].prevTime();
    const Real  prev_pres_time = state[Press_Type].prevTime();
    //
    // Compute viscous terms at level n.
    // Ensure reasonable values in 1 grow cell.  Here, do extrap for
    // c-f/phys boundary, since we have no interpolator fn, also,
    // preserve extrap for corners at periodic/non-periodic intersections.
    //
    MultiFab visc_terms(grids,nComp,1);

    if (be_cn_theta != 1.0)
    {
	getViscTerms(visc_terms,Xvel,nComp,prev_time);
    }
    else
    {
	visc_terms.setVal(0);
    }
    //
    // Set up the timestep estimation.
    //
    Real cflgrid,u_max[3];
    Real cflmax = 1.0e-10;
    comp_cfl    = (level == 0) ? cflmax : comp_cfl;

    FArrayBox tforces;

    Array<int> bndry[BL_SPACEDIM];

    MultiFab Gp(grids,BL_SPACEDIM,1);

    getGradP(Gp, prev_pres_time);
    
    FArrayBox* null_fab = 0;

    showMF("mac",Gp,"pv_Gp",level);
#if 0
    showMF("mac",get_old_data(State_Type),"pv_S_old",level);
#endif
    showMF("mac",*rho_ptime,"pv_rho_old",level);
    showMF("mac",visc_terms,"pv_visc_terms",level);

#ifndef NDEBUG
    for (int dir=0; dir<BL_SPACEDIM; ++dir) {
        u_mac[dir].setVal(0.);
    }
    MultiFab Force(grids,BL_SPACEDIM,hyp_grow);
    Force.setVal(0);
#endif

    for (FillPatchIterator U_fpi(*this,visc_terms,hyp_grow,prev_time,State_Type,Xvel,BL_SPACEDIM)
#ifdef MOREGENGETFORCE
	     , S_fpi(*this,visc_terms,1,prev_time,State_Type,Density,NUM_SCALARS);
	 S_fpi.isValid() && U_fpi.isValid();
	 ++S_fpi, ++U_fpi
#else
             ; U_fpi.isValid();
	 ++U_fpi
#endif
	)
    {
        const int i = U_fpi.index();

#ifdef GENGETFORCE
        getForce(tforces,i,1,Xvel,BL_SPACEDIM,prev_time,(*rho_ptime)[i]);
#elif MOREGENGETFORCE
	if (ParallelDescriptor::IOProcessor() && getForceVerbose)
	    std::cout << "---" << std::endl << "A - Predict velocity:" << std::endl << " Calling getForce..." << std::endl;
        getForce(tforces,i,1,Xvel,BL_SPACEDIM,prev_time,U_fpi(),S_fpi(),0);
#else
	getForce(tforces,i,1,Xvel,BL_SPACEDIM,(*rho_ptime)[i]);
#endif		 
        //
        // Test velocities, rho and cfl.
        //
        cflgrid  = godunov->test_u_rho(U_fpi(),(*rho_ptime)[i],grids[i],dx,dt,u_max);
        cflmax   = std::max(cflgrid,cflmax);
        comp_cfl = std::max(cflgrid,comp_cfl);
        //
        // Compute the total forcing.
        //
        godunov->Sum_tf_gp_visc(tforces,0,visc_terms[i],0,Gp[i],0,(*rho_ptime)[i],0);

#ifndef NDEBUG
        Force[U_fpi].copy(tforces,0,0,BL_SPACEDIM);
#endif

        D_TERM(bndry[0] = getBCArray(State_Type,i,0,1);,
               bndry[1] = getBCArray(State_Type,i,1,1);,
               bndry[2] = getBCArray(State_Type,i,2,1);)

        godunov->Setup(grids[i], dx, dt, 1,
                       D_DECL(*null_fab,*null_fab,*null_fab),
                       D_DECL(bndry[0].dataPtr(),bndry[1].dataPtr(),bndry[2].dataPtr()),
                       D_DECL(U_fpi(),U_fpi(),U_fpi()), D_DECL(0,1,2),
                       tforces,0);

        godunov->ComputeUmac(grids[i], dx, dt, 
                             u_mac[0][i], bndry[0].dataPtr(),
                             u_mac[1][i], bndry[1].dataPtr(),
#if (BL_SPACEDIM == 3)
                             u_mac[2][i], bndry[2].dataPtr(),
#endif
                             U_fpi(), tforces);
    }

    showMF("mac",u_mac[0],"pv_umac0",level);
    showMF("mac",u_mac[1],"pv_umac1",level);
#if BL_SPACEDIM==3
    showMF("mac",u_mac[2],"pv_umac2",level);
#endif
#ifndef NDEBUG
    showMF("mac",Force,"pv_force",level);
#endif

    Real tempdt = std::min(change_max,cfl/cflmax);

    ParallelDescriptor::ReduceRealMin(tempdt);

    return dt*tempdt;
}
#endif

Real
HeatTransfer::advance (Real time,
		       Real dt,
		       int  iteration,
		       int  ncycle)
{

    is_predictor = true;
    updateFluxReg = false;
    if (sdc_iterMAX == 0)
    {
      updateFluxReg = true;
    }

    if (level == 0)
    {
        crse_dt = dt;
        int thisLevelStep = parent->levelSteps(0);
        FORT_SET_COMMON(&time,&thisLevelStep);
    }

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "SDC Advancing level " << level
                  << " : starting time = " << time
                  << " with dt = "         << dt << '\n';
    }

    advance_setup(time,dt,iteration,ncycle);

    if (do_check_divudt)
        checkTimeStep(dt);
    
    MultiFab& S_new = get_new_data(State_Type);
    MultiFab& S_old = get_old_data(State_Type);

    //
    // Compute traced states for normal comp of velocity at half time level.
    //
    Real dt_test = 0.0, dummy = 0.0;    
    dt_test = predict_velocity(dt,dummy);    

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::advance(): at start of time step\n";

    temperature_stats(S_old);

    const Real prev_time = state[State_Type].prevTime();
    const Real cur_time  = state[State_Type].curTime();
    //
    // Do MAC projection and update edge velocities.
    //
    if (do_mac_proj) 
    {
        int havedivu = 1;
        create_mac_rhs(Forcing,time,dt,nGrowAdvForcing); // Use MF laying around
        mac_project(time,dt,S_old,&Forcing,havedivu);
    }

    if (do_mom_diff == 0)
        velocity_advection(dt);
    //
    // Build a copy of rho(tn) with grow cells for use in the diffusion solves
    make_rho_prev_time();

    //
    // Compute Dn and DDn (based on state at tn)
    //  (Note that coeffs at tn and tnp1 were intialized in _setup)
    //
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Computing Dn and DDn (SDC predictor) \n";

    compute_differential_diffusion_terms(Dn,DDn,prev_time,dt);

    //
    // Compute A (advection terms) with F = Dn + R
    //
    for (MFIter mfi(Forcing); mfi.isValid(); ++mfi) 
    {
        FArrayBox& f = Forcing[mfi];
        const FArrayBox& d = Dn[mfi];
        const FArrayBox& dd = DDn[mfi];
        const FArrayBox& r = get_old_data(RhoYdot_Type)[mfi];
        const Box gbox = Box(mfi.validbox()).grow(nGrowAdvForcing);
        
        f.copy(d,gbox,0,gbox,0,nspecies+1); // add Dn to RhoY and RhoH
        f.plus(r,gbox,gbox,0,0,nspecies); // add R to RhoY, no contribution for RhoH
        f.plus(dd,gbox,gbox,0,nspecies,1); // add DDn to RhoH forcing
    }
    Forcing.FillBoundary(0,nspecies+1);
    geom.FillPeriodicBoundary(Forcing,0,nspecies+1,true);

    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "A (SDC predictor)\n";

    compute_scalar_advection_fluxes_and_divergence(Forcing,dt);
    scalar_advection_update(dt, Density, RhoH);

    make_rho_curr_time();

    // 
    // Solve for Dhat, diffuse with F = A + R + 0.5*(Dn + Dhat)
    //                                        + 0.5*(DDn + DDnp1)
    // NOTE: 0.5*DDnp1 added to F within dd_update after RhoY solve but before RhoH solve
    //
    if (verbose && ParallelDescriptor::IOProcessor())
      std::cout << "Dhat (SDC predictor) \n";

    Real theta = 0.5; // c-n solve
    for (MFIter mfi(Forcing); mfi.isValid(); ++mfi) 
    {
        const Box& box = mfi.validbox();
        FArrayBox& f = Forcing[mfi];
        const FArrayBox& a = (*aofs)[mfi];
        const FArrayBox& dn = Dn[mfi];
        const FArrayBox& ddn = DDn[mfi];
        const FArrayBox& r = get_old_data(RhoYdot_Type)[mfi];
        
        f.copy(dn,box,0,box,0,nspecies+1); // copy Dn into RhoY and RhoH
        f.plus(ddn,box,box,0,nspecies,1); // add DDn to RhoH, no contribution for RhoY
        f.mult(0.5);
        f.plus(a,box,box,first_spec,0,nspecies+1); // add A into RhoY and RHoH
        f.plus(r,box,box,0,0,nspecies); // R[RhoH] == 0
    }

    // we add theta_enthalpy*DDnp1 after RhoY solve to forcing for RhoH
    Real theta_enthalpy = theta;
    
    // advection-diffusion solve with theta=0.5, theta_enthalpy=0.5
    differential_diffusion_update(Forcing,0,theta,Dhat,0,DDnp1,theta_enthalpy);

    // 
    // Compute R, react with F = A + 0.5*(Dn + Dhat)
    //                             + 0.5*(DDn + DDnp1)
    // 
    for (MFIter mfi(Forcing); mfi.isValid(); ++mfi) 
    {
        const Box& box = mfi.validbox();
        FArrayBox& f = Forcing[mfi];
        const FArrayBox& a = (*aofs)[mfi];
        const FArrayBox& dn = Dn[mfi];
        const FArrayBox& ddn = DDn[mfi];
        const FArrayBox& dhat = Dhat[mfi];
        const FArrayBox& ddnp1 = DDnp1[mfi];
        
        f.copy(dn,box,0,box,0,nspecies+1); // copy Dn into RhoY and RhoH
	f.plus(dhat,box,box,0,0,nspecies+1); // add Dhat into RhoY and RhoH
	f.plus(ddn,box,box,0,nspecies,1); // add DDn to RhoH, no contribution for RhoY
	f.plus(ddnp1,box,box,0,nspecies,1); // add DDnp1 to RhoH, no contribution for RhoY
        f.mult(0.5);
	f.plus(a,box,box,first_spec,0,nspecies+1); // add A into RhoY and RhoH
    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "R (SDC predictor) \n";

    advance_chemistry(S_old,S_new,dt,Forcing,0);

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "DONE WITH R (SDC predictor) \n";

    temperature_stats(S_new);

    for (int sdc_iter=0; sdc_iter<sdc_iterMAX; ++sdc_iter)
    {
        is_predictor = false;
        if (sdc_iter == sdc_iterMAX-1)
	{
	  updateFluxReg = true;
	}

        // Compute updated coeffs at tnp1
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "Update np1 coeffs (SDC corrector " << sdc_iter << ") \n";

        calcDiffusivity(cur_time);

	//
	// Compute Dn and DDn (based on state at tnpn)
	// iteratively lagged
	//
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "Computing Dnp1 and DDnp1 (SDC corrector " << sdc_iter << ")\n";

        compute_differential_diffusion_terms(Dnp1,DDnp1,cur_time,dt);

        //
        // Compute A (advection terms) with F = Dn + R
        //
        for (MFIter mfi(Forcing); mfi.isValid(); ++mfi) 
        {
      	    FArrayBox& f = Forcing[mfi];
            const FArrayBox& dn = Dn[mfi];
            const FArrayBox& ddn = DDn[mfi];
            const FArrayBox& r = get_new_data(RhoYdot_Type)[mfi];
	    const Box gbox = Box(mfi.validbox()).grow(nGrowAdvForcing);
            
            f.copy(dn,gbox,0,gbox,0,nspecies+1); // add Dn to RhoY and RhoH
            f.plus(r,gbox,gbox,0,0,nspecies); // add R to RhoY, no contribution for RhoH
            f.plus(ddn,gbox,gbox,0,nspecies,1); // add DDn to RhoH forcing
        }
        Forcing.FillBoundary(0,nspecies+1);
        geom.FillPeriodicBoundary(Forcing,0,nspecies,true);

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "A (SDC corrector " << sdc_iter << ")\n";

        compute_scalar_advection_fluxes_and_divergence(Forcing,dt);
        scalar_advection_update(dt, Density, RhoH);

        make_rho_curr_time();
        // 
        // Compute Dhat, diffuse with F 
	//                 = A + R + 0.5(Dn + Dnp1) - Dnp1 + Dhat + 0.5(DDn + DDnp1)
        //                 = A + R + 0.5(Dn - Dnp1) + Dhat + 0.5(DDn + DDnp1)
	// NOTE: Here we use 0.5*DDnp1 from the previous iteration
        // 
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "Dhat (SDC corrector " << sdc_iter << ")\n";

        Real sdc_theta = 1.0; // backward-Euler type solve
        for (MFIter mfi(Forcing); mfi.isValid(); ++mfi) 
        {
            const Box& box = mfi.validbox();
            FArrayBox& f = Forcing[mfi];
            const FArrayBox& a = (*aofs)[mfi];
            const FArrayBox& r = get_new_data(RhoYdot_Type)[mfi];
            const FArrayBox& dn = Dn[mfi];
            const FArrayBox& ddn = DDn[mfi];
            const FArrayBox& dnp1 = Dnp1[mfi];
	    const FArrayBox& ddnp1 = DDnp1[mfi];
            
            f.copy(dn,box,0,box,0,nspecies+1); // copy Dn into RhoY and RhoH
            f.minus(dnp1,box,box,0,0,nspecies+1); // subtract Dnp1 from RhoY and RhoH
            f.plus(ddn  ,box,box,0,nspecies,1); // add DDn to RhoH, no contribution for RhoY
	    f.plus(ddnp1,box,box,0,nspecies,1); // add DDnp1 to RhoH, no contribution for RhoY
            f.mult(0.5);
            f.plus(a,box,box,first_spec,0,nspecies+1); // add A into RhoY and RhoH
            f.plus(r,box,box,0,0,nspecies); // no reactions for RhoH
        }

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "Dhat (SDC corrector " << sdc_iter << ")\n";

	// Do not recompute enthalpy diffusion terms
        theta_enthalpy = -1;

	// advection-diffusion solve with theta=1, theta_enthalpy=-1
        differential_diffusion_update(Forcing,0,sdc_theta,Dhat,0,DDnp1,theta_enthalpy);

        // 
        // Compute R (F = A + 0.5(Dn - Dnp1 + DDn + DDnp1) + Dhat )
        // 
        for (MFIter mfi(Forcing); mfi.isValid(); ++mfi) 
        {
            const Box& box = mfi.validbox();
            FArrayBox& f = Forcing[mfi];
            const FArrayBox& a = (*aofs)[mfi];
            const FArrayBox& dn = Dn[mfi];
            const FArrayBox& dnp1 = Dnp1[mfi];
            const FArrayBox& dhat = Dhat[mfi];
            const FArrayBox& ddn = DDn[mfi];
            const FArrayBox& ddnp1 = DDnp1[mfi];

            f.copy(dn,box,0,box,0,nspecies+1); // copy Dn into RhoY and RhoH
            f.minus(dnp1,box,box,0,0,nspecies+1); // subtract Dnp1 from RhoY and RhoH
            f.plus(ddn  ,box,box,0,nspecies,1); // add DDn to RhoH, no contribution for RhoY
            f.plus(ddnp1,box,box,0,nspecies,1); // add DDnp1 to RhoH, no contribution for RhoY
            f.mult(0.5);
	    f.plus(dhat,box,box,0,0,nspecies+1); // add Dhat to RhoY and RHoH
            f.plus(a,box,box,first_spec,0,nspecies+1); // add A to RhoY and RhoH
        }

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "R (SDC corrector " << sdc_iter << ")\n";

        advance_chemistry(S_old,S_new,dt,Forcing,0);

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "DONE WITH R (SDC corrector " << sdc_iter << ")\n";

    }

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << " SDC iterations complete \n";

    calcDiffusivity(cur_time);
    calcViscosity(cur_time,dt,iteration,ncycle);
    //
    // Set the dependent value of RhoRT to be the thermodynamic pressure.  By keeping this in
    // the state, we can use the average down stuff to be sure that RhoRT_avg is avg(RhoRT),
    // not ave(Rho)avg(R)avg(T), which seems to give the p-relax stuff in the mac Rhs troubles.
    //
    setThermoPress(cur_time);

    calc_divu(time+dt, dt, get_new_data(Divu_Type));

    if (!NavierStokes::initial_step && level != parent->finestLevel())
    {
        //
        // Set new divu to old div + dt*dsdt_old where covered by fine.
        //
        BoxArray crsndgrids = getLevel(level+1).grids;

        crsndgrids.coarsen(fine_ratio);
            
        MultiFab& divu_new = get_new_data(Divu_Type);
        MultiFab& divu_old = get_old_data(Divu_Type);
        MultiFab& dsdt_old = get_old_data(Dsdt_Type);
            
        for (MFIter mfi(divu_new); mfi.isValid();++mfi)
        {
            std::vector< std::pair<int,Box> > isects = crsndgrids.intersections(mfi.validbox());

            for (int i = 0, N = isects.size(); i < N; i++)
            {
                const Box& ovlp = isects[i].second;

                divu_new[mfi].copy(dsdt_old[mfi],ovlp,0,ovlp,0,1);
                divu_new[mfi].mult(dt,ovlp,0,1);
                divu_new[mfi].plus(divu_old[mfi],ovlp,0,0,1);
            }
        }
    }
        
    calc_dsdt(time, dt, get_new_data(Dsdt_Type));

    if (NavierStokes::initial_step)
        MultiFab::Copy(get_old_data(Dsdt_Type),get_new_data(Dsdt_Type),0,0,1,0);
    //
    // Add the advective and other terms to get velocity (or momentum) at t^{n+1}.
    //
    velocity_update(dt);
    //
    // Increment rho average.
    //
    if (!initial_step)
    {
        if (level > 0)
        {
            Real alpha = 1.0/Real(ncycle);
            if (iteration == ncycle)
                alpha = 0.5/Real(ncycle);
            incrRhoAvg(alpha);
        }
        //
        // Do a level project to update the pressure and velocity fields.
        //
        level_projector(dt,time,iteration);

        if (level > 0 && iteration == 1) p_avg->setVal(0);
    }

    advance_cleanup(dt,iteration,ncycle);
    //
    // Update estimate for allowable time step.
    //
    dt_test = std::min(dt_test, estTimeStep());

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "HeatTransfer::advance(): at end of time step\n";

    temperature_stats(S_new);

    return dt_test;
}

void
HeatTransfer::create_mac_rhs (MultiFab& rhs, Real time, Real dt, int nGrow)
{
    NavierStokes::create_mac_rhs(rhs,time,dt,nGrow);
    showMF("mac",rhs,"mac_rhs0_",level);

    if (dt > 0.0) 
    {
        // BL_ASSERT(nGrow==0);  // FIXME: This does not seem to do grow cells correctly
        MultiFab  dpdt(grids,1,nGrow);
        dpdt.setVal(0.0,nGrow);
	calc_dpdt(time,dt,dpdt,u_mac);
        dpdt.FillBoundary();
        geom.FillPeriodicBoundary(dpdt,0,1);
        MultiFab::Add(rhs,dpdt,0,0,1,nGrow);
    }

    showMF("mac",rhs,"mac_rhs1_",level);
}

void
HeatTransfer::advance_chemistry (MultiFab&       mf_old,
                                 MultiFab&       mf_new,
                                 Real            dt,
                                 const MultiFab& Force,
                                 int             nCompF)
{
    const Real strt_time = ParallelDescriptor::second();

    if (hack_nochem)
    {
      BoxLib::Error("HeatTransfer.cpp: hack_nochem = T not supported yet for SDC");
    }
    else
    {
        Real p_amb, dpdt_factor;
        FORT_GETPAMB(&p_amb, &dpdt_factor);
        const Real Patm = p_amb / P1atm_MKS;

        FArrayBox* chemDiag = 0;

        if (do_not_use_funccount)
        {
	    MultiFab tmp;

            tmp.define(mf_old.boxArray(), 1, 0, mf_old.DistributionMap(), Fab_allocate);

            FArrayBox tmpf;
            for (MFIter Smfi(mf_old); Smfi.isValid(); ++Smfi)
            {
       	        FArrayBox& rYn = mf_new[Smfi];
       	        FArrayBox& rHn = mf_new[Smfi];
       	        FArrayBox& Tn  = mf_new[Smfi];
       	        const FArrayBox& rYo = mf_old[Smfi];
       	        const FArrayBox& rHo = mf_old[Smfi];
       	        const FArrayBox& To  = mf_old[Smfi];

                const Box& box = Smfi.validbox();
		FArrayBox& fc = tmp[Smfi];
		const FArrayBox& frc = Force[Smfi];
		FArrayBox& rYdot = get_new_data(RhoYdot_Type)[Smfi];

                if (plot_reactions &&
                    BoxLib::intersect(mf_old.boxArray(),auxDiag["REACTIONS"]->boxArray()).size() != 0)
                {
                    chemDiag = &( (*auxDiag["REACTIONS"])[Smfi] );
                }

                getChemSolve().solveTransient_sdc(rYn,rHn,Tn,  rYo,rHo,To,  frc, rYdot, fc, box, 
                                                  first_spec, RhoH, Temp, dt, Patm, chemDiag);
                
            }
            //
            // When ngrow>0 this does NOT properly update FuncCount_Type since parallel
            // copy()s do not touch ghost cells.  We'll ignore this since we're not using
            // the FuncCount_Type anyway.
            //
	    get_new_data(FuncCount_Type).copy(tmp);

	    //
	    // Ensure consistent grow cells
	    //    
	    MultiFab& React = get_new_data(RhoYdot_Type);
	    
	    if (React.nGrow() > 0)
	      {
		const int nc = nspecies;
		for (MFIter mfi(React); mfi.isValid(); ++mfi)
		  {
		    FArrayBox& R = React[mfi];
		    const Box& box = mfi.validbox();
		    FORT_VISCEXTRAP(R.dataPtr(),ARLIM(R.loVect()),ARLIM(R.hiVect()),
				    box.loVect(),box.hiVect(),&nc);
		  }
		React.FillBoundary(0,nspecies);
		//
		// Note: this is a special periodic fill in that we want to
		// preserve the extrapolated grow values when periodic --
		// usually we preserve only valid data.  The scheme relies on
		// the fact that there is computable data in the "non-periodic"
		// grow cells (produced via VISCEXTRAP above)
		//
		geom.FillPeriodicBoundary(React,0,nspecies,true);
	      }

        }
        else
        {
	  BoxLib::Error("HeatTransfer.cpp: do_not_use_funccount = F not supported yet for SDC");
        }
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::advance_chemistry time: " << run_time << '\n';
    }
}

void
HeatTransfer::compute_scalar_advection_fluxes_and_divergence (MultiFab& Force,
                                                              Real      dt)
{
    //
    // Compute -Div(advective fluxes)  [ which is -aofs in NS, BTW ... careful...
    //
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... computing advection terms\n";
    const Real* dx        = geom.CellSize();
    const Real prev_time = state[State_Type].prevTime();

    //
    // Gather info necesary to build transverse velocities
    // (stored internal to Godunov)
    // 
    MultiFab DivU(grids,1,nGrowAdvForcing);
    create_mac_rhs(DivU,prev_time,dt,nGrowAdvForcing);

    MultiFab Gp, VelViscTerms;
    const int use_forces_in_trans = godunov->useForcesInTrans();
    if (use_forces_in_trans || (do_mom_diff == 1))
    {
        VelViscTerms.define(grids,BL_SPACEDIM,nGrowAdvForcing,Fab_allocate);
        getViscTerms(VelViscTerms,Xvel,BL_SPACEDIM,prev_time);

        Gp.define(grids,BL_SPACEDIM,nGrowAdvForcing,Fab_allocate);
        getGradP(Gp, state[Press_Type].prevTime());
    }

    showMF("dd",DivU,"dd_divu_in_aofs",level);

    FArrayBox tforces, volume, area[BL_SPACEDIM], tvelforces, jdivu;
    const int  nState          = desc_lst[State_Type].nComp();
    for (FillPatchIterator S_fpi(*this,DivU,HYP_GROW,prev_time,State_Type,0,nState);
         S_fpi.isValid();
         ++S_fpi)
    {
        const Box& box = S_fpi.validbox();
        const int i = S_fpi.index();
        const FArrayBox& S = S_fpi();
        const FArrayBox& divu = DivU[S_fpi];
        const Box gbox = S_fpi().box();
        tforces.resize(gbox,1);

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
            geom.GetFaceArea(area[dir],grids,i,dir,GEOM_GROW);
        geom.GetVolume(volume,grids,i,GEOM_GROW);

        if (use_forces_in_trans || (do_mom_diff == 1))
        {
            NavierStokes::getForce(tvelforces,i,nGrowAdvForcing,Xvel,BL_SPACEDIM,
#ifdef GENGETFORCE
				   prev_time,
#endif		 
				   S,Density);

            godunov->Sum_tf_gp_visc(tvelforces,0,VelViscTerms[i],0,Gp[i],0,S,Density);
        }
        //
        // Set up the workspace internal to Godunov
        //
        Array<int> u_bc[BL_SPACEDIM];
        D_TERM(u_bc[0] = getBCArray(State_Type,i,0,1);,
               u_bc[1] = getBCArray(State_Type,i,1,1);,
               u_bc[2] = getBCArray(State_Type,i,2,1);)
        godunov->BuildWorkSpace(box,dx,dt);
        godunov->ComputeTransverVelocities(box,dx,dt,
                                           D_DECL(u_bc[0].dataPtr(),u_bc[1].dataPtr(),u_bc[2].dataPtr()),
                                           D_DECL(S,S,S), D_DECL(0,1,2), tvelforces, 0);
        //
        // Get edge states for Rho, RhoY, RhoH, Tracer
        //
        (*aofs)[i].setVal(0,box,Density,nspecies+3);
        for (int d=0; d<BL_SPACEDIM; ++d)
        {
            (*EdgeState[d])[i].setVal(0,(*EdgeState[d])[i].box(),Density,nspecies+3);
            (*EdgeFlux[d])[i].setVal(0,(*EdgeFlux[d])[i].box(),Density,nspecies+3);
        }        
	// loop over RhoY, RhoH, and Tracer
	// we construct density by summing RhoY
        for (int comp = 0 ; comp < nspecies+2 ; comp++)
        {
	  int state_ind = first_spec + comp;
	  Array<int> bc = getBCArray(State_Type,i,state_ind,1);
	  int iconserv = 1;

	  if (state_ind != Trac )
          {
	    godunov->edge_states_fpu(box, dx, dt,
				     D_DECL(u_mac[0][i],u_mac[1][i],u_mac[2][i]), D_DECL(0,0,0), 
				     D_DECL((*EdgeState[0])[i],(*EdgeState[1])[i],(*EdgeState[2])[i]),
				     D_DECL(state_ind,state_ind,state_ind),
				     S,state_ind,Force[S_fpi],comp,divu,0,state_ind,bc.dataPtr(),iconserv);
	  }

	  for (int d=0; d<BL_SPACEDIM; ++d)
	    (*EdgeFlux[d])[i].copy((*EdgeState[d])[i],state_ind,state_ind,1);

	  int avcomp = 0;
	  int ucomp = 0;
	  // Compute Div(flux.Area), return Area-scaled fluxes
	  godunov->ComputeAofs(grids[i],
                               D_DECL(area[0],area[1],area[2]),D_DECL(avcomp,avcomp,avcomp),
			       D_DECL(u_mac[0][i],u_mac[1][i],u_mac[2][i]),D_DECL(ucomp,ucomp,ucomp),
			       D_DECL((*EdgeFlux[0])[i],(*EdgeFlux[1])[i],(*EdgeFlux[2])[i]),
			       D_DECL(state_ind,state_ind,state_ind), volume, avcomp,
			       (*aofs)[i], state_ind, iconserv);
	  
	  // Accumulate rho flux divergence, rho on edges, and rho flux on edges
	  if (state_ind >= first_spec && state_ind < first_spec+nspecies)
          {
	    (*aofs)[i].plus((*aofs)[i],state_ind,Density,1);
	    for (int d=0; d<BL_SPACEDIM; ++d)
	      {
		(*EdgeState[d])[i].plus((*EdgeState[d])[i],state_ind,Density,1);
		(*EdgeFlux[d])[i].plus((*EdgeFlux[d])[i],state_ind,Density,1);
	      }
	  }

	  // AJN FLUXREG
	  // We have just computed the advective fluxes.  If updateFluxReg=T,
	  // ADD advective fluxes to ADVECTIVE flux register
	  if (do_reflux && updateFluxReg && level > 0)
	  {

	    // update the flux register for state_ind
	    for (int d = 0; d < BL_SPACEDIM; d++)
	      advflux_reg->FineAdd((*EdgeFlux[d])[i],d,i,state_ind,state_ind,1,dt);

	    // now that density has been constructed by summing rhoY, 
	    // update the flux register for density
	    if (state_ind == first_spec+nspecies-1)
	    {
	      for (int d = 0; d < BL_SPACEDIM; d++)
		advflux_reg->FineAdd((*EdgeFlux[d])[i],d,i,Density,Density,1,dt);
	    }

	  }
	}

	 // NOTE: Changes sense of aofs here so that d/dt ~ aofs...be sure we use our own update
	 //  function
        (*aofs)[i].mult(-1,Density,nspecies+3);
    }

    // update flux register for rho, rhoY, rhoh, and tracer
    if (do_reflux && updateFluxReg && level < parent->finestLevel())
    {
        for (int d = 0; d < BL_SPACEDIM; d++)
	{
	  getAdvFluxReg(level+1).CrseInit((*EdgeFlux[d]),d,Density,Density,nspecies+3,-dt,FluxRegister::ADD);
	}
    }
}

//
// An enum to clean up mac_sync...questionable usefullness
//
enum SYNC_SCHEME {ReAdvect, UseEdgeState, Other};

void
HeatTransfer::mac_sync ()
{
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "... mac_sync\n";

    const Real strt_time = ParallelDescriptor::second();

    int        sigma;
    const int  finest_level   = parent->finestLevel();
    const int  ngrids         = grids.size();
    const Real prev_time      = state[State_Type].prevTime();
    const Real cur_time       = state[State_Type].curTime();
    const Real prev_pres_time = state[Press_Type].prevTime();
    const Real dt             = parent->dtLevel(level);
    MultiFab*  Rh             = get_rho_half_time();
    // will hold q^{n+1,p} * (delta rho)^sync for conserved quantities
    // as defined before Eq (18) in DayBell:2000.  Note that in the paper, 
    // Eq (18) is missing Y_m^{n+1,p} * (delta rho)^sync in the RHS
    // and Eq (19) is missing the h^{n+1,p} * (delta rho)^sync in the RHS
    MultiFab*  DeltaSsync     = 0;

    //
    // Compute the corrective pressure used to compute U^{ADV,corr} in mac_sync_compute
    //
    mac_projector->mac_sync_solve(level,dt,Rh,fine_ratio);

    if (do_reflux)
    {
        MultiFab& S_new = get_new_data(State_Type);

	Array<SYNC_SCHEME> sync_scheme(NUM_STATE,ReAdvect);

        if (do_mom_diff == 1)
          for (int i=0; i<BL_SPACEDIM; ++i)
            sync_scheme[i] = UseEdgeState;

        for (int i=BL_SPACEDIM; i<NUM_STATE; ++i)
            sync_scheme[i] = UseEdgeState;
        
        Array<int> incr_sync(NUM_STATE,0);
        for (int i=0; i<sync_scheme.size(); ++i)
            if (sync_scheme[i] == ReAdvect)
                incr_sync[i] = 1;

	//
	// After solving for mac_sync_phi in mac_sync_solve(), we
	// can now do the sync advect step in mac_sync_compute().
	// This consists of two steps
	//
	// 1. compute U^{ADV,corr} as the gradient of mac_sync_phi
	// 2. add -D^MAC ( U^{ADV,corr} * rho * q)^{n+1/2} ) to flux registers
	//

	// velocities
        if (do_mom_diff == 0) 
        {
            mac_projector->mac_sync_compute(level,u_mac,Vsync,Ssync,Rh,
                                            (level > 0) ? &getAdvFluxReg(level) : 0,
                                            advectionType,prev_time,
                                            prev_pres_time,dt,NUM_STATE,
                                            be_cn_theta,
                                            modify_reflux_normal_vel,
                                            do_mom_diff,
                                            incr_sync);
        }
        else
        {
            for (int comp=0; comp<BL_SPACEDIM; ++comp)
            {
                if (sync_scheme[comp]==UseEdgeState)
                {
                    mac_projector->mac_sync_compute(level,Vsync,comp,
                                                    comp,EdgeState, comp,Rh,
                                                    (level>0 ? &getAdvFluxReg(level):0),
                                                    advectionType,modify_reflux_normal_vel,dt);
                }
            }
        }

	// scalars
        for (int comp=BL_SPACEDIM; comp<NUM_STATE; ++comp)
        {
            if (sync_scheme[comp]==UseEdgeState)
            {
                int s_ind = comp - BL_SPACEDIM;
                //
                // This routine does a sync advect step for a single 
                // scalar component. The half-time edge states are passed in.
                // This routine is useful when the edge states are computed
                // in a physics-class-specific manner. (For example, as they are
                // in the calculation of div rho U h = div U sum_l (rho Y)_l h_l(T)).
                //
                mac_projector->mac_sync_compute(level,Ssync,comp,s_ind,
                                                EdgeState,comp,Rh,
                                                (level>0 ? &getAdvFluxReg(level):0),
                                                advectionType,modify_reflux_normal_vel,dt);
            }
        }
        
        Ssync->mult(dt,Ssync->nGrow());

        sync_setup(DeltaSsync);

        //
        // For all conservative variables Q (other than density)
	// set DeltaSsync = q^{n+1,p} * (delta rho)^sync
	// subtract q^{n+1,p} * (delta rho)^sync from Ssync
        //
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const int  i   = mfi.index();
            const Box& grd = grids[i];

            int iconserved = -1;

            for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
            {
                if (istate != Density && advectionType[istate] == Conservative)
                {
                    iconserved++;
             
		    FArrayBox delta_ssync;
                    delta_ssync.resize(grd,1);
                    delta_ssync.copy(S_new[i],grd,istate,grd,0,1); // delta_ssync = (rho*q)^{n+1,p}
                    delta_ssync.divide(S_new[i],grd,Density,0,1); // delta_ssync = q^{n+1,p}
                    FArrayBox& s_sync = (*Ssync)[i]; // Ssync = RHS of Eq (18), (19) without the q^{n+1,p} * (delta rho)^sync terms
                    delta_ssync.mult(s_sync,grd,Density-BL_SPACEDIM,0,1); // delta_ssync = q^{n+1,p} * (delta rho)^sync
                    (*DeltaSsync)[i].copy(delta_ssync,grd,0,grd,iconserved,1); // DeltaSsync = q^{n+1,p} * (delta rho)^sync
                    s_sync.minus(delta_ssync,grd,0,istate-BL_SPACEDIM,1); // Ssync = Ssync - q^{n+1,p} * (delta rho)^sync
                }
            }
        }
        //
        // Now, increment density.
        //
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();
            S_new[i].plus((*Ssync)[i],grids[i],Density-BL_SPACEDIM,Density,1);
        }
        make_rho_curr_time();

        const int numscal = NUM_STATE - BL_SPACEDIM;
        //
        // Set do_diffuse_sync to 0 for debugging reasons only.
        //
        if (do_mom_diff == 1)
        {
            for (MFIter Vsyncmfi(*Vsync); Vsyncmfi.isValid(); ++Vsyncmfi)
            {
                const int i    = Vsyncmfi.index();
                const Box vbox = (*rho_ctime).box(i);

                D_TERM((*Vsync)[i].divide((*rho_ctime)[i],vbox,0,Xvel,1);,
                       (*Vsync)[i].divide((*rho_ctime)[i],vbox,0,Yvel,1);,
                       (*Vsync)[i].divide((*rho_ctime)[i],vbox,0,Zvel,1););
            }
        }

        if (do_diffuse_sync)
        {
	    MultiFab** beta;
	    diffusion->allocFluxBoxesLevel(beta);
            if (is_diffusive[Xvel])
            {
                int rho_flag = (do_mom_diff == 0) ? 1 : 3;
                getViscosity(beta, cur_time);
                diffusion->diffuse_Vsync(Vsync,dt,be_cn_theta,Rh,rho_flag,beta);
            }
	    
	    if (!unity_Le 
		&& nspecies>0 
		&& do_add_nonunityLe_corr_to_rhoh_adv_flux 
		&& !do_mcdd)
	    {
		//
		// Diffuse the species syncs such that sum(SpecDiffSyncFluxes) = 0
	        // After exiting, SpecDiffusionFluxnp1 should contain rhoD grad (delta Y)^sync
		// Also, Ssync for species should contain rho^{n+1} * (delta Y)^sync
		differential_spec_diffuse_sync(dt);

                MultiFab Soln(grids,1,1);
                const Real cur_time  = state[State_Type].curTime();
                const Real a = 1.0;     // Passed around, but not used
                Real rhsscale;          //  -ditto-
                const int rho_flag = 2; // FIXME: Messy assumption
                MultiFab *alpha=0;      //  -ditto-
                MultiFab **fluxSC, **fluxNULN, **rhoh_visc;
                diffusion->allocFluxBoxesLevel(fluxSC,0,1);
                diffusion->allocFluxBoxesLevel(fluxNULN,0,nspecies);
                diffusion->allocFluxBoxesLevel(rhoh_visc,0,1);

                const int nGrow    = 1; // Size to grow fil-patched fab for T below
                const int dataComp = 0; // coeffs loaded into 0-comp for all species
                  
		// get lambda/cp
                getDiffusivity(rhoh_visc, cur_time, RhoH, 0, 1);

		// compute lambda/cp grad (delta Y_m^sync)
                for (int comp = 0; comp < nspecies; ++comp)
                {
                    const Real b     = be_cn_theta;
                    const int  sigma = first_spec + comp;
                    //
                    //  start by getting lambda/cp.Grad(delta Y^sync)
                    //   (note: neg of usual diff flux)
                    //
                    ViscBndry      visc_bndry;
                    ABecLaplacian* visc_op;

                    visc_op = diffusion->getViscOp(sigma,a,b,cur_time,
                                                   visc_bndry,Rh,
                                                   rho_flag,&rhsscale,rhoh_visc,dataComp,
                                                   alpha,dataComp);

                    visc_op->maxOrder(diffusion->maxOrder());

		    // copy rho^{n+1} * (delta Y)^sync into Soln
                    MultiFab::Copy(Soln,*Ssync,sigma-BL_SPACEDIM,0,1,0);

		    // divide Soln by rho
                    for (MFIter Smfi(Soln); Smfi.isValid(); ++Smfi)
		    {
		      Soln[Smfi].divide(S_new[Smfi],Smfi.validbox(),Density,0,1);
		    }

		    // compute lambda/cp.Grad(delta Y^sync) and weight
		    visc_op->compFlux(D_DECL(*fluxSC[0],*fluxSC[1],*fluxSC[2]),Soln);
                    for (int d = 0; d < BL_SPACEDIM; ++d)
                        fluxSC[d]->mult(-b/geom.CellSize()[d]);
                    //
                    // Here, get fluxNULN = (lambda/cp - rho.D)Grad(delta Ysync)
                    //                    = lambda/cp.Grad(delta Ysync) + SpecSyncDiffFlux
                    //
                    for (int d = 0; d < BL_SPACEDIM; ++d)
                    {
                        MFIter SDF_mfi(*SpecDiffusionFluxnp1[d]);

                        for ( ; SDF_mfi.isValid(); ++SDF_mfi)
                        {
                            FArrayBox& fluxSC_fab   = (*fluxSC[d])[SDF_mfi];
                            FArrayBox& fluxNULN_fab = (*fluxNULN[d])[SDF_mfi];
                            FArrayBox& SDF_fab = (*SpecDiffusionFluxnp1[d])[SDF_mfi];
                            const Box& ebox    = SDF_mfi.validbox();
			    // copy in (delta Gamma)^sync
                            fluxNULN_fab.copy(SDF_fab,ebox,comp,ebox,comp,1);
			    // add in (lambda/cp) grad (delta Y^sync)
                            fluxNULN_fab.plus(fluxSC_fab,ebox,0,comp,1);
                        }
                    }
                    delete visc_op;
                }

                Soln.clear();
                diffusion->removeFluxBoxesLevel(fluxSC);
                diffusion->removeFluxBoxesLevel(rhoh_visc);
                //
                // Multiply fluxi by h_i (let FLXDIV routine below sum up the fluxes)
                //
                FArrayBox eTemp, h, volume;

                for (FillPatchIterator Tnew_fpi(*this,S_new,nGrow,cur_time,State_Type,Temp,1);
                     Tnew_fpi.isValid();
                     ++Tnew_fpi)
                {
                    const int i    = Tnew_fpi.index();
                    const Box& box = Tnew_fpi.validbox();

                    for (int d = 0; d < BL_SPACEDIM; ++d)
                    {
                        const Box ebox = BoxLib::surroundingNodes(box,d);
                        eTemp.resize(ebox,1);
                        FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).lo(d));
                        FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Temp).hi(d));
                        
                        center_to_edge_fancy(Tnew_fpi(),eTemp,BoxLib::grow(box,BoxLib::BASISV(d)),
                                             0,0,1,geom.Domain(),bc_lo,bc_hi);
                        
                        h.resize(ebox,nspecies);
                        getChemSolve().getHGivenT(h,eTemp,ebox,0,0);

			// multiply fluxNULN by h_m
                        (*fluxNULN[d])[i].mult(h,ebox,0,0,nspecies);
                    }
                }

		// add the NULN fluxes to the RHS of the (delta h)^sync diffusion solve
		// afterwards, the entire RHS should be ready
                for (MFIter Ssync_mfi(*Ssync); Ssync_mfi.isValid(); ++Ssync_mfi)
                {
                    const int i      = Ssync_mfi.index();
                    FArrayBox& syncn = (*Ssync)[i];
                    const FArrayBox& synco = (*Ssync)[i];
                    const Box& box = Ssync_mfi.validbox();

                    geom.GetVolume(volume,grids,i,GEOM_GROW);
                    //
                    // Multiply by dt*dt, one to make it extensive, and one because
                    // Ssync multiplied above by dt, need same units here.
                    //
                    const Real mult = dt*dt;
                    const int sigmaRhoH = RhoH - BL_SPACEDIM; // RhoH comp in Ssync
		    
                    FORT_INCRWEXTFLXDIV(box.loVect(), box.hiVect(),
                                        (*fluxNULN[0])[i].dataPtr(),
                                        ARLIM((*fluxNULN[0])[i].loVect()),
                                        ARLIM((*fluxNULN[0])[i].hiVect()),
                                        (*fluxNULN[1])[i].dataPtr(),
                                        ARLIM((*fluxNULN[1])[i].loVect()),
                                        ARLIM((*fluxNULN[1])[i].hiVect()),
#if BL_SPACEDIM == 3
                                        (*fluxNULN[2])[i].dataPtr(),
                                        ARLIM((*fluxNULN[2])[i].loVect()),
                                        ARLIM((*fluxNULN[2])[i].hiVect()),
#endif
                                        synco.dataPtr(sigmaRhoH),
                                        ARLIM(synco.loVect()),
                                        ARLIM(synco.hiVect()),
                                        syncn.dataPtr(sigmaRhoH),
                                        ARLIM(syncn.loVect()),
                                        ARLIM(syncn.hiVect()),
                                        volume.dataPtr(),
                                        ARLIM(volume.loVect()),
                                        ARLIM(volume.hiVect()),
                                        &nspecies, &mult);
                }

                diffusion->removeFluxBoxesLevel(fluxNULN);
            }
            else if (nspecies>0 && do_mcdd)
            {
                mcdd_diffuse_sync(dt);
            }

	    MultiFab **flux;
            diffusion->allocFluxBoxesLevel(flux);

            for (sigma = 0; sigma < numscal; sigma++)
            {
                int rho_flag = 0;
                int do_viscsyncflux = do_reflux;
		const int state_ind = BL_SPACEDIM + sigma;
		//
		// To diffuse, or not?
		// (1) Density, no
		// (2) RhoH...if diffusive
		// (3) Trac...if diffusive
		// (4) Spec:
		//    (a) if Le==1, and spec diffusive
		//    (b) if Le!=1, do differential diffusion instead, done above
		// (5) Temp, no (set instead by RhoH to Temp)
                //
		const bool is_spec = state_ind<=last_spec && state_ind>=first_spec;
                int do_it
		    =  state_ind!=Density 
		    && state_ind!=Temp
		    && is_diffusive[state_ind]
		    && !(is_spec && !unity_Le)
		    && !(do_mcdd && (is_spec || state_ind==RhoH));
		
		if (do_it && (is_spec || state_ind==RhoH))
		    rho_flag = 2;

                if (do_it)
                {
                    MultiFab* alpha = 0;
                    getDiffusivity(beta, cur_time, state_ind, 0, 1);
                    
		    // on entry, Ssync = RHS for (delta h)^sync diffusive solve
		    // on exit, Ssync = rho^{n+1} * (delta h)^sync
		    // on exit, flux = coeff * grad phi
		    diffusion->diffuse_Ssync(Ssync,sigma,dt,be_cn_theta,Rh,
					     rho_flag,flux,0,beta,0,alpha);
		    if (do_viscsyncflux && level > 0)
		    {
			for (MFIter mfi(*Ssync); mfi.isValid(); ++mfi)
			{
			    const int i=mfi.index();
			    for (int d=0; d<BL_SPACEDIM; ++d)
                                getViscFluxReg().FineAdd((*flux[d])[i],d,i,0,state_ind,1,dt);
			}
		    }
                }
            }
            diffusion->removeFluxBoxesLevel(flux);
	    diffusion->removeFluxBoxesLevel(beta);
        }
        //
        // For all conservative variables Q (other than density)
        // increment sync by (sync_for_rho)*q_presync.
	// Before this loop, Ssync holds rho^{n+1} (delta phi)^sync
	// DeltaSsync holds (delta rho)^sync phi^p
        //
        for (MFIter mfi(*Ssync); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            int iconserved = -1;

            for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
            {
                if (istate != Density && advectionType[istate] == Conservative)
                {
                    iconserved++;

                    (*Ssync)[i].plus((*DeltaSsync)[i],grids[i],iconserved,istate-BL_SPACEDIM,1);
                }
            }
        }
        sync_cleanup(DeltaSsync);
        //
        // Increment the state (for all but rho, since that was done above)
        //
        for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            for (int sigma = 0; sigma < numscal; sigma++)
            {
                if (!(BL_SPACEDIM+sigma == Density))
                {
                    S_new[i].plus((*Ssync)[i],grids[i],sigma,BL_SPACEDIM+sigma,1);
                }
            }
        }
        //
        // Recompute temperature and rho R T after the mac_sync.
        //
        RhoH_to_Temp(S_new);
        setThermoPress(cur_time);
        //
        // Get boundary conditions.
        //
        Real mult = 1.0;
        Array<int*>         sync_bc(grids.size());
        Array< Array<int> > sync_bc_array(grids.size());
        for (int i = 0; i < ngrids; i++)
        {
            sync_bc_array[i] = getBCArray(State_Type,i,Density,numscal);
            sync_bc[i]       = sync_bc_array[i].dataPtr();
        }
        //
        // Interpolate the sync correction to the finer levels.
        //
        IntVect ratio = IntVect::TheUnitVector();
        for (int lev = level+1; lev <= finest_level; lev++)
        {
            ratio                   *= parent->refRatio(lev-1);
            HeatTransfer& fine_level = getLevel(lev);
            MultiFab& S_new_lev      = fine_level.get_new_data(State_Type);
            //
            // New way of interpolating syncs to make sure mass is conserved
            // and to ensure freestream preservation for species & temperature.
            //
            const BoxArray& fine_grids = S_new_lev.boxArray();
            const int nghost           = S_new_lev.nGrow();
            MultiFab increment(fine_grids, numscal, nghost);
            increment.setVal(0,nghost);
            //
            // Note: we use the lincc_interp (which_interp==3) for density,
            // rho*h and rho*Y, cell_cons_interp for everything else. Doing
            // so is needed for freestream preservation of Y and T.  The setting
            // which_interp=5 calls the lincc_interp, but then follows with a
            // "protection" step to be sure that all the components but the
            // first and the last have their sync adjusted to try to preserve
            // positivity after the sync is applied.  The density sync is then
            // adjusted to be the sum of the species syncs.
            //
            // HACK note: Presently, the species mass syncs are redistributed 
            //            without consequence to the enthalpy sync.  Clearly
            //            the species carry enthalphy, so the enthalph sync should
            //            be adjusted as well.  Note yet sure how to do this correctly.
            //            Punt for now...
            //
            const SyncInterpType which_interp = CellConsProt_T;

            const int nComp = 2+nspecies;

            SyncInterp(*Ssync, level, increment, lev, ratio, 
                       Density-BL_SPACEDIM, Density-BL_SPACEDIM, nComp, 1, mult, 
                       sync_bc.dataPtr(), which_interp, Density);

            if (have_trac)
                SyncInterp(*Ssync, level, increment, lev, ratio, 
                           Trac-BL_SPACEDIM, Trac-BL_SPACEDIM, 1, 1, mult, 
                           sync_bc.dataPtr());

            if (have_rhort)
                SyncInterp(*Ssync, level, increment, lev, ratio, 
                           RhoRT-BL_SPACEDIM, RhoRT-BL_SPACEDIM, 1, 1, mult, 
                           sync_bc.dataPtr());

            SyncInterp(*Ssync, level, increment, lev, ratio, 
                       Temp-BL_SPACEDIM, Temp-BL_SPACEDIM, 1, 1, mult, 
                       sync_bc.dataPtr());

            if (do_set_rho_to_species_sum)
            {
                increment.setVal(0,Density-BL_SPACEDIM,1,0);

                for (int istate = first_spec; istate <= last_spec; istate++)
                { 
                    for (MFIter mfi(increment); mfi.isValid(); ++mfi)
                    {
                        int i = mfi.index();
                        increment[i].plus(increment[i],fine_grids[i],
					  istate-BL_SPACEDIM,Density-BL_SPACEDIM,1);
                    }
                }
            }

            for (MFIter mfi(increment); mfi.isValid(); ++mfi)
            {
                int i = mfi.index();
                S_new_lev[i].plus(increment[i],fine_grids[i],0,Density,numscal);
            }
            fine_level.make_rho_curr_time();
            fine_level.incrRhoAvg(increment,Density-BL_SPACEDIM,1.0);
            //
            // Recompute temperature and rho R T after interpolation of the mac_sync correction
            //   of the individual quantities rho, Y, T.
            //
            RhoH_to_Temp(S_new_lev);
            fine_level.setThermoPress(cur_time);
        }
        //
        // Average down Trac = rho R T after interpolation of the mac_sync correction
        //   of the individual quantities rho, Y, T.
        //
        for (int lev = finest_level-1; lev >= level; lev--)
        {
            HeatTransfer&   fine_lev = getLevel(lev+1);
            const BoxArray& fgrids   = fine_lev.grids;
            
            HeatTransfer&   crse_lev = getLevel(lev);
            const BoxArray& cgrids   = crse_lev.grids;
            const IntVect&  fratio   = crse_lev.fine_ratio;
            
            MultiFab& S_crse = crse_lev.get_new_data(State_Type);
            MultiFab& S_fine = fine_lev.get_new_data(State_Type);

            const int pComp = (have_rhort ? RhoRT : Trac);
            crse_lev.NavierStokes::avgDown(cgrids,fgrids,S_crse,S_fine,
                                           lev,lev+1,pComp,1,fratio);
        }
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer:mac_sync(): lev: " << level << ", time: " << run_time << '\n';
    }
}

void
HeatTransfer::mcdd_diffuse_sync(Real dt)
{
    //
    // Compute the increment such that the composite equations are satisfied.
    // In the level advance, we solved
    //  (rho.Y)^* - (rho.Y)^n + theta.dt.L(Y)^* + (1-theta).dt.L(Y)^n = -dt.Div(rho.U.Y)
    // Here, we solve
    //  (rho.Y)^np1 - (rho.Y)^* + theta.dt(L(Y)^np1 - L(Y)^*) = -dt.Div(corr) = Ssync
    // where corr is the area.time weighted correction to add to the coarse to get fine
    //
    if (hack_nospecdiff || hack_nomcddsync)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! skipping mcdd sync diffusion " << '\n';

        for (int d=0; d<BL_SPACEDIM; ++d)
            SpecDiffusionFluxnp1[d]->setVal(0);

        return;            
    }
    //
    // TODO - merge in new CrseInit() calls.
    //
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "...doing mcdd sync diffusion..." << '\n';

    const Real time = state[State_Type].curTime();
    const int nGrow = 0;
    MultiFab Rhs(grids,nspecies+1,nGrow,Fab_allocate);
    const int dCompY = 0;
    const int dCompH = nspecies;
    PArray<MultiFab> fluxpre(BL_SPACEDIM,PArrayManage);
    for (int dir=0; dir<BL_SPACEDIM; ++dir)
        fluxpre.set(dir,new MultiFab(BoxArray(grids).surroundingNodes(dir),
                                     nspecies+1,0));
    
    // Rhs = Ssync + dt*theta*L(phi_presync) + (rhoY)_presync
    compute_mcdd_visc_terms(Rhs,time,nGrow,DDOp::DD_RhoH,&fluxpre);

    const int spec_Ssync_sComp = first_spec - BL_SPACEDIM;
    const int rhoh_Ssync_sComp = RhoH - BL_SPACEDIM;
    MultiFab& S_new = get_new_data(State_Type);
    
    for (MFIter mfi(Rhs); mfi.isValid(); ++mfi)
    {
        FArrayBox& rhs = Rhs[mfi];
        const FArrayBox& snew = S_new[mfi];
        const FArrayBox& sync = (*Ssync)[mfi];
        const Box& box = mfi.validbox();

        rhs.mult(be_cn_theta*dt,box,0,nspecies+1);

        rhs.plus(sync,box,spec_Ssync_sComp,dCompY,nspecies);
        rhs.plus(sync,box,rhoh_Ssync_sComp,dCompH,1);

        rhs.plus(snew,box,first_spec,dCompY,nspecies);
        rhs.plus(snew,box,RhoH,dCompH,1);
    }

    // Save a copy of the pre-sync state
    MultiFab::Copy(*Ssync,S_new,first_spec,spec_Ssync_sComp,nspecies,0);
    MultiFab::Copy(*Ssync,S_new,RhoH,rhoh_Ssync_sComp,1,0);

    PArray<MultiFab> fluxpost(BL_SPACEDIM,PArrayManage);
    for (int dir=0; dir<BL_SPACEDIM; ++dir)
        fluxpost.set(dir,new MultiFab(BoxArray(grids).surroundingNodes(dir),
                                      nspecies+1,nGrow)); //stack H on Y's
    
    // Solve the system
    //mcdd_solve(Rhs,dCompY,dCompH,fluxpost,dCompY,dCompH,dt);
    BoxLib::Error("Fix the sync diffuse!");

    // Set Ssync to hold the resulting increment
    for (MFIter mfi(S_new); mfi.isValid(); ++mfi)
    {
        FArrayBox& sync = (*Ssync)[mfi];
        const FArrayBox& snew = S_new[mfi];
        const Box& box = mfi.validbox();
        
        sync.negate(box,spec_Ssync_sComp,nspecies);
        sync.negate(box,rhoh_Ssync_sComp,1);
        
        sync.plus(snew,box,first_spec,spec_Ssync_sComp,nspecies);
        sync.plus(snew,box,RhoH,rhoh_Ssync_sComp,1);
    }

    if (do_reflux)
    {
        {
            FArrayBox fluxinc;
            for (int d = 0; d < BL_SPACEDIM; d++)
            {
                for (MFIter mfi(fluxpost[d]); mfi.isValid(); ++mfi)
                {
                    const FArrayBox& fpost = fluxpost[d][mfi];
                    const FArrayBox& fpre = fluxpre[d][mfi];
                    const Box& ebox = fpost.box();

                    fluxinc.resize(ebox,nspecies+1);
                    fluxinc.copy(fpost,ebox,0,ebox,0,nspecies+1);
                    fluxinc.minus(fpre,ebox,0,0,nspecies+1);
                    fluxinc.mult(be_cn_theta);
                    if (level < parent->finestLevel())
                    {
                        //
                        // Note: The following inits do not trash each other
                        // since the component ranges don't overlap...
                        //
                        getLevel(level+1).getViscFluxReg().CrseInit(fluxinc,ebox,
                                                                    d,dCompY,first_spec,
                                                                    nspecies,-dt);
                        getLevel(level+1).getViscFluxReg().CrseInit(fluxinc,ebox,
                                                                    d,dCompH,RhoH,
                                                                    1,-dt);
                    }
                    if (level > 0)
                    {
                        getViscFluxReg().FineAdd(fluxinc,d,mfi.index(),
                                                 dCompY,first_spec,nspecies,dt);
                        getViscFluxReg().FineAdd(fluxinc,d,mfi.index(),
                                                 dCompH,RhoH,1,dt);
                    }
                }
            }
        }
        if (level < parent->finestLevel())
            getLevel(level+1).getViscFluxReg().CrseInitFinish();
    }
    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "...finished differential sync diffusion..." << '\n';    
}

void
HeatTransfer::differential_spec_diffuse_sync (Real dt)
{
  
  // Diffuse the species syncs such that sum(SpecDiffSyncFluxes) = 0
  // After exiting, SpecDiffusionFluxnp1 should contain rhoD grad (delta Y)^sync
  // Also, Ssync for species should contain rho^{n+1} * (delta Y)^sync

    if (hack_nospecdiff)
    {
        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "... HACK!!! skipping spec sync diffusion " << '\n';

        for (int d=0; d<BL_SPACEDIM; ++d)
            SpecDiffusionFluxnp1[d]->setVal(0.0,0,nspecies);

        return;            
    }

    const Real strt_time = ParallelDescriptor::second();

    if (verbose && ParallelDescriptor::IOProcessor())
        std::cout << "Doing differential sync diffusion ..." << '\n';
    //
    // Do implicit c-n solve for each scalar...but dont reflux.
    // Save the fluxes, coeffs and source term, we need 'em for later
    // Actually, since Ssync multiplied by dt in mac_sync before
    // a call to this routine (I hope...), Ssync is in units of s,
    // not the "usual" ds/dt...lets convert it (divide by dt) so
    // we can use a generic flux adjustment function
    //
    const Real cur_time = state[State_Type].curTime();
    MultiFab **betanp1;
    diffusion->allocFluxBoxesLevel(betanp1,0,nspecies);
    getDiffusivity(betanp1, cur_time, first_spec, 0, nspecies); // rhoD

    MultiFab Rhs(grids,nspecies,0);
    const int spec_Ssync_sComp = first_spec - BL_SPACEDIM;

    // Rhs and Ssync contain the RHS of DayBell:2000 Eq (18)
    // with the additional -Y_m^{n+1,p} * (delta rho)^sync term
    // Copy this into Rhs; we will need this later since we overwrite SSync in the solves
    MultiFab::Copy(Rhs,*Ssync,spec_Ssync_sComp,0,nspecies,0);

    //
    // Some standard settings
    //
    const Array<int> rho_flag(nspecies,2);
    const MultiFab* alpha = 0;
    MultiFab** fluxSC;
    diffusion->allocFluxBoxesLevel(fluxSC,0,1);

    const MultiFab* RhoHalftime = get_rho_half_time();

    for (int sigma = 0; sigma < nspecies; ++sigma)
    {
        //
        // Here, we use Ssync as a source in units of s, as expected by diffuse_Ssync
        // (i.e., ds/dt ~ d(Ssync)/dt, vs. ds/dt ~ Rhs in diffuse_scalar).  This was
        // apparently done to mimic diffuse_Vsync, which does the same, because the
        // diffused result is an acceleration, not a velocity, req'd by the projection.
        //
	const int ssync_ind = first_spec + sigma - Density;
                    
	// on entry, Ssync = RHS for (delta Ytilde)^sync diffusive solve
	// on exit, Ssync = rho^{n+1} * (delta Ytilde)^sync
	// on exit, fluxSC = rhoD grad (delta Ytilde)^sync
	diffusion->diffuse_Ssync(Ssync,ssync_ind,dt,be_cn_theta,
				 RhoHalftime,rho_flag[sigma],fluxSC,0,
                                 betanp1,sigma,alpha);
	//
	// Pull fluxes into flux array
	// this is the rhoD grad (delta Ytilde)^sync terms in DayBell:2000 Eq (18)
	//
	for (int d=0; d<BL_SPACEDIM; ++d)
	{
	  MultiFab::Copy(*SpecDiffusionFluxnp1[d],*fluxSC[d],0,sigma,1,0);
	}
    }
    diffusion->removeFluxBoxesLevel(fluxSC);
    //
    // Modify update/fluxes to preserve flux sum = 0
    // (Be sure to pass the "normal" looking Rhs to this generic function)
    //

    // need to correct SpecDiffusionFluxnp1 to contain rhoD grad (delta Y)^sync
    bool grow_cells_already_filled = false;
    adjust_spec_diffusion_fluxes(cur_time,grow_cells_already_filled);

    // need to correct Ssync to contain rho^{n+1} * (delta Y)^sync
    // Do this by setting
    // Ssync = "RHS from diffusion solve" + (dt/2)*div(delta Gamma)
    //
    // Recompute update with adjusted diffusion fluxes
    // 

    FArrayBox update, volume, efab[BL_SPACEDIM];

    for (MFIter mfi(*Ssync); mfi.isValid(); ++mfi)
    {
	int        iGrid = mfi.index();
	const Box& box   = mfi.validbox();

	// copy corrected (delta gamma) on edges into efab
        for (int d=0; d<BL_SPACEDIM; ++d)
        {
            const Box ebox = BoxLib::surroundingNodes(box,d);
            efab[d].resize(ebox,nspecies);
            
            efab[d].copy((*SpecDiffusionFluxnp1[d])[mfi],ebox,0,ebox,0,nspecies);
            efab[d].mult(be_cn_theta);
        }

        update.resize(box,nspecies);
	update.setVal(0.);
        geom.GetVolume(volume,grids,iGrid,GEOM_GROW);

	// is this right? - I'm not sure what the scaling on SpecDiffusionFluxnp1 is
        Real scale = -dt;

	// take divergence of (delta gamma) and put it in update
	// update = scale * div(flux) / vol
	// we want update to contain (dt/2) div (delta gamma)
	FORT_FLUXDIV(box.loVect(), box.hiVect(),
                     update.dataPtr(),  
		     ARLIM(update.loVect()),  ARLIM(update.hiVect()),
                     efab[0].dataPtr(), ARLIM(efab[0].loVect()), ARLIM(efab[0].hiVect()),
                     efab[1].dataPtr(), ARLIM(efab[1].loVect()), ARLIM(efab[1].hiVect()),
#if BL_SPACEDIM == 3
                     efab[2].dataPtr(), ARLIM(efab[2].loVect()), ARLIM(efab[2].hiVect()),
#endif
                     volume.dataPtr(),  ARLIM(volume.loVect()),  ARLIM(volume.hiVect()),
                     &nspecies,&scale);

	// add RHS from diffusion solve
	update.plus(Rhs[iGrid],box,0,0,nspecies);

	// Ssync = "RHS from diffusion solve" - dt/2) div (delta gamma)
	(*Ssync)[mfi].copy(update,box,0,box,first_spec-BL_SPACEDIM,nspecies);

    }
    diffusion->removeFluxBoxesLevel(betanp1);

    Rhs.clear();
    //
    // Do refluxing AFTER flux adjustment
    //
    if (do_reflux && level > 0)
    {
      for (int d=0; d<BL_SPACEDIM; ++d)
      {
	for (MFIter fmfi(*SpecDiffusionFluxnp1[d]); fmfi.isValid(); ++fmfi)
	{
	  const int i=fmfi.index();
	  getViscFluxReg().FineAdd((*SpecDiffusionFluxnp1[d])[i],d,i,0,first_spec,nspecies,dt);
	}
      }
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::differential_spec_diffuse_sync(): time: " << run_time << '\n';
    }
}

void
HeatTransfer::reflux ()
{
    // no need to reflux if this is the finest level
    if (level == parent->finestLevel()) return;

    const Real strt_time = ParallelDescriptor::second();

    BL_ASSERT(do_reflux);

    //
    // First do refluxing step.
    //
    FluxRegister& fr_adv  = getAdvFluxReg(level+1);
    FluxRegister& fr_visc = getViscFluxReg(level+1);
    const Real    dt_crse = parent->dtLevel(level);
    const Real    scale   = 1.0/dt_crse;
    //
    // It is important, for do_mom_diff == 0, to do the viscous
    //   refluxing first, since this will be divided by rho_half
    //   before the advective refluxing is added.  In the case of
    //   do_mom_diff == 1, both components of the refluxing will
    //   be divided by rho^(n+1) in NavierStokes::level_sync.
    //
    MultiFab volume;

    geom.GetVolume(volume,grids,GEOM_GROW);

    // take divergence of diffusive flux registers into cell-centered RHS
    fr_visc.Reflux(*Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    if (true)
        fr_visc.Reflux(*Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    const MultiFab* RhoHalftime = get_rho_half_time();

    if (do_mom_diff == 0) 
    {
        for (MFIter mfi(*Vsync); mfi.isValid(); ++mfi)
        {
            const int i = mfi.index();

            D_TERM((*Vsync)[i].divide((*RhoHalftime)[i],grids[i],0,Xvel,1);,
                   (*Vsync)[i].divide((*RhoHalftime)[i],grids[i],0,Yvel,1);,
                   (*Vsync)[i].divide((*RhoHalftime)[i],grids[i],0,Zvel,1););
        }
    }

    FArrayBox tmp;

    // for any variables that used non-conservative advective differencing,
    // divide the sync by rhohalf
    for (MFIter mfi(*Ssync); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();

        tmp.resize(grids[i],1);
        tmp.copy((*RhoHalftime)[i],0,0,1);
        tmp.invert(1);

        for (int istate = BL_SPACEDIM; istate < NUM_STATE; istate++)
        {
            if (advectionType[istate] == NonConservative)
            {
                const int sigma = istate -  BL_SPACEDIM;

                (*Ssync)[i].mult(tmp,0,sigma,1);
            }
        }
    }

    // take divergence of advective flux registers into cell-centered RHS
    fr_adv.Reflux(*Vsync,volume,scale,0,0,BL_SPACEDIM,geom);
    fr_adv.Reflux(*Ssync,volume,scale,BL_SPACEDIM,0,NUM_STATE-BL_SPACEDIM,geom);

    BoxArray baf = getLevel(level+1).boxArray();

    baf.coarsen(fine_ratio);
    //
    // This is necessary in order to zero out the contribution to any
    // coarse grid cells which underlie fine grid cells.
    //
    std::vector< std::pair<int,Box> > isects;

    isects.reserve(27);

    for (MFIter mfi(*Vsync); mfi.isValid(); ++mfi)
    {
        baf.intersections(grids[mfi.index()],isects);

        for (int i = 0, N = isects.size(); i < N; i++)
        {
            (*Vsync)[mfi].setVal(0,isects[i].second,0,BL_SPACEDIM);
            (*Ssync)[mfi].setVal(0,isects[i].second,0,NUM_STATE-BL_SPACEDIM);
        }
    }

    if (verbose)
    {
        const int IOProc   = ParallelDescriptor::IOProcessorNumber();
        Real      run_time = ParallelDescriptor::second() - strt_time;

        ParallelDescriptor::ReduceRealMax(run_time,IOProc);

        if (ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::Reflux(): time: " << run_time << '\n';
    }
}

void
HeatTransfer::calcViscosity (const Real time,
			     const Real dt,
			     const int  iteration,
			     const int  ncycle)
{
    const TimeLevel whichTime = which_time(State_Type, time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    compute_vel_visc(time, whichTime == AmrOldTime ? viscn_cc : viscnp1_cc);
}

void
HeatTransfer::calcDiffusivity (const Real time)
{
    calcDiffusivity(time,false);
}

void
HeatTransfer::calcDiffusivity (const Real time,
                               bool       do_VelVisc)
{
    if (do_mcdd) return;

    const TimeLevel whichTime = which_time(State_Type, time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    const int  nGrow           = 1;
    const int  offset          = BL_SPACEDIM + 1; // No diffusion coeff for vels or rho
    MultiFab&  diff            = (whichTime == AmrOldTime) ? (*diffn_cc) : (*diffnp1_cc);
    FArrayBox tmp, bcen;
    for (FillPatchIterator Rho_and_spec_fpi(*this,diff,nGrow,time,State_Type,Density,nspecies+1),
             Temp_fpi(*this,diff,nGrow,time,State_Type,Temp,1);
         Rho_and_spec_fpi.isValid() && Temp_fpi.isValid();
         ++Rho_and_spec_fpi, ++Temp_fpi)
    {
        FArrayBox& Tfab = Temp_fpi();
        FArrayBox& RYfab = Rho_and_spec_fpi();
         const Box& gbox = RYfab.box();
        //
        // Convert from RhoY_l to Y_l
        //
        tmp.resize(gbox,1);
        tmp.copy(RYfab,0,0,1);
        tmp.invert(1);

	//	for (int n = 1; n < nspecies+1; n++)
	//	  RYfab.mult(tmp,0,n,1);

        const int  vflag   = do_VelVisc;
        const int nc_bcen = nspecies+2; // rhoD + lambda + mu
        int       dotemp  = 1;
        bcen.resize(gbox,nc_bcen);
        
        FORT_SPECTEMPVISC(gbox.loVect(),gbox.hiVect(),
                          ARLIM(Tfab.loVect()),ARLIM(Tfab.hiVect()),Tfab.dataPtr(),
                          ARLIM(RYfab.loVect()),ARLIM(RYfab.hiVect()),RYfab.dataPtr(1),
                          ARLIM(bcen.loVect()),ARLIM(bcen.hiVect()),bcen.dataPtr(),
                          &nc_bcen, &P1atm_MKS, &dotemp, &vflag);
        
        FArrayBox& Dfab = diff[Rho_and_spec_fpi];
        Dfab.copy(bcen,0,first_spec-offset,nspecies);
        Dfab.copy(bcen,nspecies,Temp-offset,1);

        if (do_VelVisc)
        {
            MultiFab& visc = (whichTime==AmrOldTime) ? (*viscn_cc) : (*viscnp1_cc);
            visc[Rho_and_spec_fpi].copy(bcen,nspecies+1,0,1);
        }


	for (int n = 1; n < nspecies+1; n++)
	  RYfab.mult(tmp,0,n,1);

        for (int icomp = RhoH; icomp <= NUM_STATE; icomp++)
        {
            if (icomp == RhoH)
            {
                const int sCompT = 0, sCompY = 1, sCompCp = RhoH-offset;
                getChemSolve().getCpmixGivenTY(Dfab,Tfab,RYfab,gbox,sCompT,sCompY,sCompCp);
                Dfab.invert(1,sCompCp,1);
                Dfab.mult(Dfab,Temp-offset,RhoH-offset,1);
            }
            else if (icomp == Trac || icomp == RhoRT)
            {
                Dfab.setVal(trac_diff_coef, gbox, icomp-offset, 1);
            }
        }
    }
    showMFsub("1D",diff,stripBox,"1D_calcD_visc",level);
}

void
HeatTransfer::getViscosity (MultiFab*  beta[BL_SPACEDIM],
                            const Real time)
{
    const TimeLevel whichTime = which_time(State_Type, time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    MultiFab* visc = (whichTime == AmrOldTime) ? viscn_cc : viscnp1_cc;

    for (MFIter viscMfi(*visc); viscMfi.isValid(); ++viscMfi)
    {
        const int i = viscMfi.index();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Density).lo(dir));
            FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(Density).hi(dir));

            center_to_edge_fancy((*visc)[viscMfi],(*beta[dir])[i],
                                 BoxLib::grow(grids[i],BoxLib::BASISV(dir)), 0, 0, 1,
                                 geom.Domain(), bc_lo, bc_hi);
        }
    }
}


void
HeatTransfer::getDiffusivity (MultiFab*  beta[BL_SPACEDIM],
                              const Real time,
                              const int  state_comp,
                              const int  dst_comp,
                              const int  ncomp)
{
    BL_ASSERT(state_comp > Density);

    const TimeLevel whichTime = which_time(State_Type, time);

    BL_ASSERT(whichTime == AmrOldTime || whichTime == AmrNewTime);

    MultiFab* diff      = (whichTime == AmrOldTime) ? diffn_cc : diffnp1_cc;
    const int offset    = BL_SPACEDIM + 1; // No diffusion coeff for vels or rho
    int       diff_comp = state_comp - offset;

    for (MFIter diffMfi(*diff); diffMfi.isValid(); ++diffMfi)
    {
        const int i = diffMfi.index();

        for (int dir = 0; dir < BL_SPACEDIM; dir++)
        {
            FPLoc bc_lo = fpi_phys_loc(get_desc_lst()[State_Type].getBC(state_comp).lo(dir));
            FPLoc bc_hi = fpi_phys_loc(get_desc_lst()[State_Type].getBC(state_comp).hi(dir));

            center_to_edge_fancy((*diff)[diffMfi],(*beta[dir])[i],
                                 BoxLib::grow(grids[i],BoxLib::BASISV(dir)), diff_comp, 
                                 dst_comp, ncomp, geom.Domain(), bc_lo, bc_hi);
        }
    }

    if (zeroBndryVisc > 0)
        zeroBoundaryVisc(beta,time,state_comp,dst_comp,ncomp);
}

void
HeatTransfer::zeroBoundaryVisc (MultiFab*  beta[BL_SPACEDIM],
                                const Real time,
                                const int  state_comp,
                                const int  dst_comp,
                                const int  ncomp) const
{
    BL_ASSERT(state_comp > Density);

    const int isrz = (int) geom.IsRZ();
    for (int dir = 0; dir < BL_SPACEDIM; dir++)
    {
        Box edom = BoxLib::surroundingNodes(geom.Domain(),dir);
        
        for (MFIter mfi(*(beta[dir])); mfi.isValid(); ++mfi)
        {
            FArrayBox& beta_fab = (*(beta[dir]))[mfi];
            const Box ebox      = BoxLib::surroundingNodes(mfi.validbox(),dir);
            FORT_ZEROVISC(beta_fab.dataPtr(dst_comp),
                          ARLIM(beta_fab.loVect()), ARLIM(beta_fab.hiVect()),
                          ebox.loVect(),  ebox.hiVect(),
                          edom.loVect(),  edom.hiVect(),
                          geom.CellSize(), geom.ProbLo(), phys_bc.vect(),
                          &dir, &isrz, &state_comp, &ncomp);
        }
    }
}

void
HeatTransfer::compute_vel_visc (Real      time,
                                MultiFab* beta)
{
    const int nGrow = beta->nGrow();

    BL_ASSERT(nGrow == 1);
    BL_ASSERT(first_spec == Density+1);

    FArrayBox tmp;

    MultiFab dummy(grids,1,0,Fab_noallocate);

    for (FillPatchIterator Temp_fpi(*this,dummy,nGrow,time,State_Type,Temp,1),
             Rho_and_spec_fpi(*this,dummy,nGrow,time,State_Type,Density,nspecies+1);
         Rho_and_spec_fpi.isValid() && Temp_fpi.isValid();
         ++Rho_and_spec_fpi, ++Temp_fpi)
    {
        const int  i            = Rho_and_spec_fpi.index();
        const Box& box          = Rho_and_spec_fpi().box();
        FArrayBox& temp         = Temp_fpi();
        FArrayBox& rho_and_spec = Rho_and_spec_fpi();
        //
        // Convert from RhoY_l to Y_l
        //
        tmp.resize(box,1);
        tmp.copy(rho_and_spec,0,0,1);
        tmp.invert(1);

        for (int n = 1; n < nspecies+1; ++n)
            rho_and_spec.mult(tmp,0,n,1);

        FORT_VELVISC(box.loVect(),box.hiVect(),
                     ARLIM(temp.loVect()),ARLIM(temp.hiVect()),temp.dataPtr(),
                     ARLIM(rho_and_spec.loVect()),ARLIM(rho_and_spec.hiVect()),
                     rho_and_spec.dataPtr(1),
                     ARLIM(tmp.loVect()),ARLIM(tmp.hiVect()),tmp.dataPtr());

        (*beta)[i].copy(tmp,0,0,1);
    }
}

void
HeatTransfer::calc_divu (Real      time,
                         Real      dt,
                         MultiFab& divu)
{
    const int nGrow = 0;
    MultiFab mcViscTerms;
    int vtCompY;
    int vtCompT;
    if (do_mcdd)
    {
        vtCompT = nspecies;
        vtCompY = 0;
        mcViscTerms.define(grids,nspecies+1,nGrow,Fab_allocate);
        compute_mcdd_visc_terms(mcViscTerms,time,nGrow,DDOp::DD_Temp);
    }
    else
    {
        vtCompT = nspecies + 1;
        vtCompY = 0;
        mcViscTerms.define(grids,nspecies+2,nGrow,Fab_allocate); // Can probably use DDnp1 here
        compute_differential_diffusion_terms(mcViscTerms,divu,time,dt); // divu has DD, not needed here
    }

    MultiFab& S = get_data(State_Type,time);
    MultiFab RhoYdot(grids,nspecies,0);
    
    if (dt > 0)
      {
	compute_instantaneous_reaction_rates(RhoYdot,S,nGrow);
      }
    else
      {
	RhoYdot.setVal(0);
      }
    
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        const Box& box = mfi.validbox();            
        FArrayBox& du = divu[mfi];
        const FArrayBox& vtY = mcViscTerms[mfi];
        const FArrayBox& vtT = mcViscTerms[mfi];
        const FArrayBox& rhoY = S[mfi];
        const FArrayBox& rYdot = RhoYdot[mfi];
        const FArrayBox& T = S[mfi];
        FORT_CALCDIVU(box.loVect(), box.hiVect(),
                      du.dataPtr(),             ARLIM(du.loVect()),    ARLIM(du.hiVect()),
                      rYdot.dataPtr(),          ARLIM(rYdot.loVect()), ARLIM(rYdot.hiVect()),
                      vtY.dataPtr(vtCompY),     ARLIM(vtY.loVect()),   ARLIM(vtY.hiVect()),
                      vtT.dataPtr(vtCompT),     ARLIM(vtT.loVect()),   ARLIM(vtT.hiVect()),
                      rhoY.dataPtr(first_spec), ARLIM(rhoY.loVect()),  ARLIM(rhoY.hiVect()),
                      T.dataPtr(Temp),          ARLIM(T.loVect()),     ARLIM(T.hiVect()));
    }
}

//
// Compute the Eulerian Dp/Dt for use in pressure relaxation.
//
void
HeatTransfer::calc_dpdt (Real      time,
                         Real      dt_,
                         MultiFab& dpdt,
                         MultiFab* u_mac)
{
    Real dt = crse_dt;

    Real p_amb, dpdt_factor;
    
    FORT_GETPAMB(&p_amb, &dpdt_factor);
    
    if (dt <= 0.0 || dpdt_factor <= 0)
    {
        dpdt.setVal(0);
        return;
    }
    
    MultiFab& S_old = get_old_data(State_Type);
    
    const int pComp = (have_rhort ? RhoRT : Trac);
    const int nGrow = 1;
    
    int sCompR, sCompT, sCompY, sCompCp, sCompMw;
    sCompR=sCompT=sCompY=sCompCp=sCompMw=0;
    
    MultiFab temp(grids,1,nGrow);
    MultiFab rhoRT(grids,1,nGrow);
    
    const MultiFab& rho = get_rho(time);

    BL_ASSERT(rho.boxArray()  == S_old.boxArray());
    BL_ASSERT(temp.boxArray() == S_old.boxArray());
    BL_ASSERT(dpdt.boxArray() == S_old.boxArray());
    
    for (FillPatchIterator Temp_fpi(*this,temp,nGrow,time,State_Type,Temp,1);
         Temp_fpi.isValid();
         ++Temp_fpi)
    {
        temp[Temp_fpi.index()].copy(Temp_fpi(),0,sCompT,1);
    }

    for (FillPatchIterator rhoRT_fpi(*this,rhoRT,nGrow,time,State_Type,pComp,1);
         rhoRT_fpi.isValid();
         ++rhoRT_fpi)
    {
        rhoRT[rhoRT_fpi.index()].copy(rhoRT_fpi(),0,0,1);
    }
    //
    // Note that state contains rho*species, so divide species by rho.
    //
    MultiFab species(grids,nspecies,nGrow);

    FArrayBox tmp;

    for (FillPatchIterator Spec_fpi(*this,species,nGrow,time,State_Type,first_spec,nspecies);
         Spec_fpi.isValid();
         ++Spec_fpi)
    {
        const int i  = Spec_fpi.index();
        const Box bx = BoxLib::grow(grids[i],nGrow);

        species[i].copy(Spec_fpi(),0,sCompY,nspecies);

        tmp.resize(bx,1);
        tmp.copy(rho[i],sCompR,0,1);
        tmp.invert(1);

        for (int ispecies = 0; ispecies < nspecies; ispecies++)
            species[i].mult(tmp,0,ispecies,1);
    }

    tmp.clear();
    
    MultiFab mwmix(grids,1,nGrow), cp(grids,1,nGrow);
    
    for (MFIter Rho_mfi(rho); Rho_mfi.isValid(); ++Rho_mfi)
    {
        const int idx = Rho_mfi.index();
        const Box box = BoxLib::grow(Rho_mfi.validbox(),nGrow);
        
        BL_ASSERT(Rho_mfi.validbox() == grids[idx]);
        
        getChemSolve().getMwmixGivenY(mwmix[idx],species[idx],box,sCompY,sCompMw);
        getChemSolve().getCpmixGivenTY(cp[idx],temp[idx],species[idx],box,sCompT,sCompY,sCompCp);
    }

    temp.clear();
    species.clear();
    //
    // Now do pressure relaxation.
    //
    const Real* dx = geom.CellSize();
    //
    //  Get gamma = c_p/(c_p-R)
    //
    MultiFab gamma(grids,1,nGrow);

    FArrayBox rmix;

    for (MFIter mfi(gamma); mfi.isValid(); ++mfi)
    {
        const int idx = mfi.index();
        gamma[idx].copy(cp[idx]);
        rmix.resize(BoxLib::grow(grids[idx],nGrow),1);
        rmix.setVal(rgas);
        rmix.divide(mwmix[idx]);
        gamma[idx].minus(rmix);
        gamma[idx].divide(cp[idx]);
        gamma[idx].invert(1.0);
    }

    rmix.clear();
    cp.clear();
    mwmix.clear();

    FArrayBox ugradp, p_denom;
    
    for (MFIter mfi(dpdt); mfi.isValid(); ++mfi)
    {
        const int i = mfi.index();
        
        dpdt[i].copy(S_old[i],grids[i],pComp,grids[i],0,1);
        dpdt[i].plus(-p_amb,grids[i]);
        dpdt[i].mult(1.0/dt,grids[i]);

        ugradp.resize(grids[i],1);
        
        const int* lo = grids[i].loVect();
        const int* hi = grids[i].hiVect();

        FORT_COMPUTE_UGRADP(rhoRT[i].dataPtr(0), 
                            ARLIM(rhoRT[i].box().loVect()), 
                            ARLIM(rhoRT[i].box().hiVect()),
                            ugradp.dataPtr(), 
                            ARLIM(ugradp.box().loVect()), 
                            ARLIM(ugradp.box().hiVect()),
                            u_mac[0][i].dataPtr(),
                            ARLIM(u_mac[0][i].box().loVect()),
                            ARLIM(u_mac[0][i].box().hiVect()),
                            u_mac[1][i].dataPtr(),
                            ARLIM(u_mac[1][i].box().loVect()),
                            ARLIM(u_mac[1][i].box().hiVect()),
#if (BL_SPACEDIM == 3)
                            u_mac[2][i].dataPtr(),
                            ARLIM(u_mac[2][i].box().loVect()),
                            ARLIM(u_mac[2][i].box().hiVect()),
#endif
                            lo,hi,dx);
        //
        // Note that this is minus because the term we want to add to S is
        //  - Dp/Dt.   We already have (p-p0)/dt = -dp/dt, so now
        // we subtract ugradp to form -Dp/dt.  
        // (Note the dp/dt term is negative because
        //  p is the current value, p0 is the value we're trying to get to,
        //  so dp/dt = (p0 - p)/dt.)
        //
        dpdt[i].minus(ugradp,grids[i],0,0,1);
        //
        // Make sure to divide by gamma *after* subtracting ugradp.
        //
        dpdt[i].divide(gamma[i],grids[i],0,0,1);

        if (dpdt_option == 0)
        {
	  dpdt[i].divide(S_old[i],grids[i],pComp,0,1);
        }
        else if (dpdt_option == 1)
        {
	  dpdt[i].divide(p_amb,grids[i],0,1);
        }
        else
        {
            const int ncomp = 1;
            p_denom.resize(grids[i],ncomp);
            p_denom.copy(S_old[i],grids[i],pComp,grids[i],0,ncomp);
            Real num_norm = dpdt[i].norm(0,0,1);
            FabMinMax(p_denom,grids[i],2*num_norm*Real_MIN,p_amb,0,ncomp);
            dpdt[i].divide(p_denom,grids[i],0,0,1);
        }
        
        dpdt[i].mult(dpdt_factor,grids[i]);
    }

    if (dpdt.nGrow() > 0)
      {
	const int nc = 1;
	for (MFIter mfi(dpdt); mfi.isValid(); ++mfi)
	  {
	    FArrayBox& dpdtfab = dpdt[mfi];
	    const Box& box = mfi.validbox();
	    FORT_VISCEXTRAP(dpdtfab.dataPtr(),ARLIM(dpdtfab.loVect()),ARLIM(dpdtfab.hiVect()),
			    box.loVect(),box.hiVect(),&nc);
	  }
	dpdt.FillBoundary(0,1);
	geom.FillPeriodicBoundary(dpdt,0,1,true);
      }
}

//
// Function to use if Divu_Type and Dsdt_Type are in the state.
//

void
HeatTransfer::calc_dsdt (Real      time,
                         Real      dt,
                         MultiFab& dsdt)
{
    MultiFab& Divu_new = get_new_data(Divu_Type);
    MultiFab& Divu_old = get_old_data(Divu_Type);

    dsdt.setVal(0);

    for (MFIter mfi(dsdt); mfi.isValid(); ++mfi)
    {
        const int k = mfi.index();
        dsdt[k].copy(Divu_new[k],grids[k]);
        dsdt[k].minus(Divu_old[k],grids[k],0,0,1);
        dsdt[k].divide(dt,grids[k],0,1);
    }
}

int
HeatTransfer::RhoH_to_Temp (FArrayBox& S,
                            const Box& box,
                            int        dominmax)
{
    BL_ASSERT(S.box().contains(box));
    //
    // Convert rho to 1/rho, rho*h to h and rho*Y to Y for this operation.
    //
    S.invert(1,box,Density,1);    
    S.mult(S,Density,RhoH,1);
    for (int spec = first_spec; spec <= last_spec; spec++)
        S.mult(S,Density,spec,1);

    int iters = HY_to_Temp_DoIt(S,S,S,box,RhoH,first_spec,Temp,htt_hmixTYP,getChemSolve());

    if (dominmax)
        FabMinMax(S, box, htt_tempmin, htt_tempmax, Temp, 1);
    //
    // Convert back to rho, rho*h and rho*Y
    //
    S.invert(1,box,Density,1);    
    S.mult(S,box,Density,RhoH,1);
    for (int spec = first_spec; spec <= last_spec; spec++)
        S.mult(S,box,Density,spec,1);

    return iters;
}

void
HeatTransfer::RhoH_to_Temp (MultiFab& S,
                            int       nGrow,
                            int       dominmax)
{
    int max_iters = 0;
    for (MFIter mfi(S); mfi.isValid(); ++mfi)
    {
        Box box = Box(mfi.validbox()).grow(nGrow);
        max_iters = std::max(max_iters, RhoH_to_Temp(S[mfi],box,dominmax));
    }

    if (verbose)
    {
        const int IOProc = ParallelDescriptor::IOProcessorNumber();

        ParallelDescriptor::ReduceIntMax(max_iters,IOProc);

        if (verbose && ParallelDescriptor::IOProcessor())
            std::cout << "HeatTransfer::RhoH_to_Temp: max_iters = " << max_iters << '\n';
    }
}

void
HeatTransfer::compute_rhohmix (Real      time,
                               MultiFab& rhohmix)
{
    const int ngrow  = 0; // We only do this on the valid region
    const int sComp  = std::min(std::min((int)Density,(int)Temp),first_spec);
    const int eComp  = std::max(std::max((int)Density,(int)Temp),first_spec+nspecies-1);
    const int nComp  = eComp - sComp + 1;
    const int sCompR = Density - sComp;
    const int sCompT = Temp - sComp;
    const int sCompY = first_spec - sComp;
    const int sCompH = 0;

    FArrayBox tmp;

    for (FillPatchIterator fpi(*this,rhohmix,ngrow,time,State_Type,sComp,nComp);
         fpi.isValid();
         ++fpi)
    {
        //
        // Convert rho*Y to Y for this operation
        //
	const int  iGrid    = fpi.index();
	const Box& validBox = grids[iGrid];
	FArrayBox& state    = fpi();

        tmp.resize(validBox,1);
        tmp.copy(state,sCompR,0,1);
        tmp.invert(1);

        for (int k = 0; k < nspecies; k++)
            state.mult(tmp,0,sCompY+k,1);
	
	getChemSolve().getHmixGivenTY(rhohmix[iGrid],state,state,validBox,
				      sCompT,sCompY,sCompH);
        //
        // Convert hmix to rho*hmix
        //
        rhohmix[iGrid].mult(state,validBox,sCompR,sCompH,1);
    }
}

void
HeatTransfer::setPlotVariables ()
{
    AmrLevel::setPlotVariables();

    const Array<std::string>& names = getChemSolve().speciesNames();

    ParmParse pp("ht");

    bool plot_ydot,plot_rhoY,plot_massFrac,plot_moleFrac,plot_conc;
    plot_ydot=plot_rhoY=plot_massFrac=plot_moleFrac=plot_conc = false;

    if (pp.query("plot_massfrac",plot_massFrac))
    {
        if (plot_massFrac)
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "Y("+names[i]+")";
                parent->addDerivePlotVar(name);
            }
        }
        else
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "Y("+names[i]+")";
                parent->deleteDerivePlotVar(name);
            }
        }
    }

    if (pp.query("plot_molefrac",plot_moleFrac))
    {
        if (plot_moleFrac)
            parent->addDerivePlotVar("molefrac");
        else
            parent->deleteDerivePlotVar("molefrac");
    }

    if (pp.query("plot_concentration",plot_conc))
    {
        if (plot_conc)
            parent->addDerivePlotVar("concentration");
        else
            parent->deleteDerivePlotVar("concentration");
    }

    if (pp.query("plot_rhoY",plot_rhoY))
    {
        if (plot_rhoY)
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "rho.Y("+names[i]+")";
                parent->addStatePlotVar(name);
            }
        }
        else
        {
            for (int i = 0; i < names.size(); i++)
            {
                const std::string name = "rho.Y("+names[i]+")";
                parent->deleteStatePlotVar(name);
            }
        }
    }

    if (verbose && ParallelDescriptor::IOProcessor())
    {
        std::cout << "\nState Plot Vars: ";

        std::list<std::string>::const_iterator li = parent->statePlotVars().begin(), end = parent->statePlotVars().end();

        for ( ; li != end; ++li)
            std::cout << *li << ' ';
        std::cout << '\n';

        std::cout << "\nDerive Plot Vars: ";

        li  = parent->derivePlotVars().begin();
        end = parent->derivePlotVars().end();

        for ( ; li != end; ++li)
            std::cout << *li << ' ';
        std::cout << '\n';
    }
}

void
HeatTransfer::writePlotFile (const std::string& dir,
                             std::ostream&  os,
                             VisMF::How     how)
{
    if ( ! Amr::Plot_Files_Output() ) return;
    //
    // Note that this is really the same as its NavierStokes counterpart,
    // but in order to add diagnostic MultiFabs into the plotfile, code had
    // to be interspersed within this function.
    //
    int i, n;
    //
    // The list of indices of State to write to plotfile.
    // first component of pair is state_type,
    // second component of pair is component # within the state_type
    //
    std::vector<std::pair<int,int> > plot_var_map;
    for (int typ = 0; typ < desc_lst.size(); typ++)
        for (int comp = 0; comp < desc_lst[typ].nComp();comp++)
            if (parent->isStatePlotVar(desc_lst[typ].name(comp)) &&
                desc_lst[typ].getType() == IndexType::TheCellType())
                plot_var_map.push_back(std::pair<int,int>(typ,comp));

    int num_derive = 0;
    std::list<std::string> derive_names;
    const std::list<DeriveRec>& dlist = derive_lst.dlist();

    for (std::list<DeriveRec>::const_iterator it = dlist.begin(), end = dlist.end();
         it != end;
         ++it)
    {
        if (parent->isDerivePlotVar(it->name()))
	{
            derive_names.push_back(it->name());
            num_derive += it->numDerive();
	}
    }

    int num_auxDiag = 0;
    for (std::map<std::string,MultiFab*>::const_iterator it = auxDiag.begin(), end = auxDiag.end();
         it != end;
         ++it)
    {
        num_auxDiag += it->second->nComp();
    }

    int n_data_items = plot_var_map.size() + num_derive + num_auxDiag;
    Real cur_time = state[State_Type].curTime();

    if (level == 0 && ParallelDescriptor::IOProcessor())
    {
        //
        // The first thing we write out is the plotfile type.
        //
        os << thePlotFileType() << '\n';

        if (n_data_items == 0)
            BoxLib::Error("Must specify at least one valid data item to plot");

        os << n_data_items << '\n';

	//
	// Names of variables -- first state, then derived
	//
	for (i =0; i < plot_var_map.size(); i++)
        {
	    int typ  = plot_var_map[i].first;
	    int comp = plot_var_map[i].second;
	    os << desc_lst[typ].name(comp) << '\n';
        }

	for (std::list<std::string>::const_iterator it = derive_names.begin(), end = derive_names.end();
             it != end;
             ++it)
        {
	    const DeriveRec* rec = derive_lst.get(*it);
	    for (i = 0; i < rec->numDerive(); i++)
                os << rec->variableName(i) << '\n';
        }
        //
        // Hack in additional diagnostics.
        //
        for (std::map<std::string,Array<std::string> >::const_iterator it = auxDiag_names.begin(), end = auxDiag_names.end();
             it != end;
             ++it)
        {
            for (i=0; i<it->second.size(); ++i)
                os << it->second[i] << '\n';
        }

        os << BL_SPACEDIM << '\n';
        os << parent->cumTime() << '\n';
        int f_lev = parent->finestLevel();
        os << f_lev << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbLo(i) << ' ';
        os << '\n';
        for (i = 0; i < BL_SPACEDIM; i++)
            os << Geometry::ProbHi(i) << ' ';
        os << '\n';
        for (i = 0; i < f_lev; i++)
            os << parent->refRatio(i)[0] << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->Geom(i).Domain() << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
            os << parent->levelSteps(i) << ' ';
        os << '\n';
        for (i = 0; i <= f_lev; i++)
        {
            for (int k = 0; k < BL_SPACEDIM; k++)
                os << parent->Geom(i).CellSize()[k] << ' ';
            os << '\n';
        }
        os << (int) Geometry::Coord() << '\n';
        os << "0\n"; // Write bndry data.
    }
    // Build the directory to hold the MultiFab at this level.
    // The name is relative to the directory containing the Header file.
    //
    static const std::string BaseName = "/Cell";

    std::string Level = BoxLib::Concatenate("Level_", level, 1);
    //
    // Now for the full pathname of that directory.
    //
    std::string FullPath = dir;
    if (!FullPath.empty() && FullPath[FullPath.length()-1] != '/')
        FullPath += '/';
    FullPath += Level;
    //
    // Only the I/O processor makes the directory if it doesn't already exist.
    //
    if (ParallelDescriptor::IOProcessor())
        if (!BoxLib::UtilCreateDirectory(FullPath, 0755))
            BoxLib::CreateDirectoryFailed(FullPath);
    //
    // Force other processors to wait till directory is built.
    //
    ParallelDescriptor::Barrier();

    if (ParallelDescriptor::IOProcessor())
    {
        os << level << ' ' << grids.size() << ' ' << cur_time << '\n';
        os << parent->levelSteps(level) << '\n';

        for (i = 0; i < grids.size(); ++i)
        {
            RealBox gridloc = RealBox(grids[i],geom.CellSize(),geom.ProbLo());
            for (n = 0; n < BL_SPACEDIM; n++)
                os << gridloc.lo(n) << ' ' << gridloc.hi(n) << '\n';
        }
        //
        // The full relative pathname of the MultiFabs at this level.
        // The name is relative to the Header file containing this name.
        // It's the name that gets written into the Header.
        //
        if (n_data_items > 0)
        {
            std::string PathNameInHeader = Level;
            PathNameInHeader += BaseName;
            os << PathNameInHeader << '\n';
        }
    }
    //
    // We combine all of the multifabs -- state, derived, etc -- into one
    // multifab -- plotMF.
    // NOTE: we are assuming that each state variable has one component,
    // but a derived variable is allowed to have multiple components.
    int       cnt   = 0;
    int       ncomp = 1;
    const int nGrow = 0;
    MultiFab  plotMF(grids,n_data_items,nGrow);
    MultiFab* this_dat = 0;
    //
    // Cull data from state variables -- use no ghost cells.
    //
    for (i = 0; i < plot_var_map.size(); i++)
    {
	int typ  = plot_var_map[i].first;
	int comp = plot_var_map[i].second;
	this_dat = &state[typ].newData();
	MultiFab::Copy(plotMF,*this_dat,comp,cnt,ncomp,nGrow);
	cnt+= ncomp;
    }
    //
    // Cull data from derived variables.
    // 
    Real plot_time;

    if (derive_names.size() > 0)
    {
	for (std::list<std::string>::const_iterator it = derive_names.begin(), end = derive_names.end();
             it != end;
             ++it) 
	{
            if (*it == "avg_pressure" ||
                *it == "gradpx"       ||
                *it == "gradpy"       ||
                *it == "gradpz") 
            {
                if (state[Press_Type].descriptor()->timeType() == 
                    StateDescriptor::Interval) 
                {
                    plot_time = cur_time;
                }
                else
                {
                    int f_lev = parent->finestLevel();
                    plot_time = getLevel(f_lev).state[Press_Type].curTime();
                }
            }
            else
            {
                plot_time = cur_time;
            } 
	    const DeriveRec* rec = derive_lst.get(*it);
	    ncomp = rec->numDerive();
	    MultiFab* derive_dat = derive(*it,plot_time,nGrow);
	    MultiFab::Copy(plotMF,*derive_dat,0,cnt,ncomp,nGrow);
	    delete derive_dat;
	    cnt += ncomp;
	}
    }
    //
    // Cull data from diagnostic multifabs.
    //
    for (std::map<std::string,MultiFab*>::const_iterator it = auxDiag.begin(), end = auxDiag.end();
         it != end;
         ++it)
    {
        int nComp = it->second->nComp();
        MultiFab::Copy(plotMF,*it->second,0,cnt,nComp,nGrow);
        cnt += nComp;
    }
    //
    // Use the Full pathname when naming the MultiFab.
    //
    std::string TheFullPath = FullPath;
    TheFullPath += BaseName;
    VisMF::Write(plotMF,TheFullPath,how);
}

MultiFab*
HeatTransfer::derive (const std::string& name,
                      Real               time,
                      int                ngrow)
{
#ifdef PARTICLES
    return ParticleDerive(name,time,ngrow);
#else
    return AmrLevel::derive(name,time,ngrow);
#endif
}

void
HeatTransfer::derive (const std::string& name,
                      Real               time,
                      MultiFab&          mf,
                      int                dcomp)
{
    AmrLevel::derive(name,time,mf,dcomp);
}

#ifdef PARTICLES
MultiFab*
HeatTransfer::ParticleDerive(const std::string& name,
                             Real               time,
                             int                ngrow)
{
    if (HTPC && name == "particle_count")
    {
        MultiFab* derive_dat = new MultiFab(grids,1,0);
        MultiFab    temp_dat(grids,1,0);
        temp_dat.setVal(0);
        HTPC->Increment(temp_dat,level);
        MultiFab::Copy(*derive_dat,temp_dat,0,0,1,0);
        return derive_dat;
    }
    else if (HTPC && name == "total_particle_count")
    {
        //
        // We want the total particle count at this level or higher.
        //
        MultiFab* derive_dat = ParticleDerive("particle_count",time,ngrow);

        IntVect trr(D_DECL(1,1,1));

        for (int lev = level+1; lev <= parent->finestLevel(); lev++)
        {
            BoxArray ba = parent->boxArray(lev);

            MultiFab temp_dat(ba,1,0);

            trr *= parent->refRatio(lev-1);

            ba.coarsen(trr);

            MultiFab ctemp_dat(ba,1,0);

            temp_dat.setVal(0);
            ctemp_dat.setVal(0);

            HTPC->Increment(temp_dat,lev);

            for (MFIter mfi(temp_dat); mfi.isValid(); ++mfi)
            {
                const FArrayBox& ffab =  temp_dat[mfi];
                FArrayBox&       cfab = ctemp_dat[mfi];
                const Box&       fbx  = ffab.box();

                BL_ASSERT(cfab.box() == BoxLib::coarsen(fbx,trr));

                for (IntVect p = fbx.smallEnd(); p <= fbx.bigEnd(); fbx.next(p))
                {
                    const Real val = ffab(p);
                    if (val > 0)
                        cfab(BoxLib::coarsen(p,trr)) += val;
                }
            }

            temp_dat.clear();

            MultiFab dat(grids,1,0);
            dat.setVal(0);
            dat.copy(ctemp_dat);

            MultiFab::Add(*derive_dat,dat,0,0,1,0);
        }

        return derive_dat;
    }
    else
    {
        return AmrLevel::derive(name,time,ngrow);
    }
}
#endif

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// XrdCpWorkLst                                                         //
//                                                                      //
// Author: Fabrizio Furano (INFN Padova, 2004)                          //
//                                                                      //
// A class implementing a list of cp to do for XrdCp                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

//   $Id$

#include <sys/types.h>
#include <dirent.h>
#include "XrdClient/XrdClientAdmin.hh"
#include "XrdClient/XrdClient.hh"


class XrdCpWorkLst {

   vecString fWorkList;
   int fWorkIt;

   XrdClientAdmin *xrda_src, *xrda_dst;

   XrdClientString fSrc, fDest;
   bool fDestIsDir, fSrcIsDir;

 public:
   
   XrdCpWorkLst();
   ~XrdCpWorkLst();

   // Sets the source path for the file copy
   int SetSrc(XrdClient **srccli, XrdClientString url,
	      XrdClientString urlopaquedata, bool do_recurse);

   // Sets the destination of the file copy
   int SetDest(XrdClient **xrddest, const char *url,
	       const char *urlopaquedata,
	       kXR_unt16 xrdopenflags);

   inline void GetDest(XrdClientString &dest, bool& isdir) {
      dest = fDest;
      isdir = fDestIsDir;
   }

   // Actually builds the worklist
   int BuildWorkList_xrd(XrdClientString url, XrdClientString opaquedata);
   int BuildWorkList_loc(DIR *dir, XrdClientString pat);

   bool GetCpJob(XrdClientString &src, XrdClientString &dest);
   
};

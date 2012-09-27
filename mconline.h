/* $Author: kobics $   */
/* $Date: 2003/06/26 07:21:18 $     */
/* $Source: /cs/phd/kobics/.CVSROOT/code/multiClass/mconline.h,v $   */
/* $Name:  $     */
/* $Locker:  $   */
/* $Revision: 6.5 $ */
/* $State: Exp $    */



#ifndef __MCONLINE_H
#define __MCONLINE_H

#include "kernel.h"
#include "redopt.h"
#include "mucutils.h"

enum MCOnlineUpdateType 
  { UPDATE_UNIFORM=0, 
    UPDATE_MAX,        /* 1 */
    UPDATE_PROP,       /* 2 */
    UPDATE_MIRA,       /* 3 */
    UPDATE_PERCEPTRON, /* 4 */ 
    UPDATE_OVR_MIRA,   /* 5 */ 
    UPDATE_MIRA1,       /* 6 not defined */ 
    UPDATE_RAND        /* 7 */ 
  };

enum MCOnlineFindType
  { FIND_MAXIMAL_MARGIN=0,
    FIND_MIMIMAL_WEIGHT,
    FIND_MAXIMAL_MARGIN_WO_EXAMPLE, /* CleanCache 1 */
    FIND_MAXIMAL_DIFF_NORM, /* CleanCache 2 */
    FIND_MAXIMAL_DIFF_NORM_NORMALIZED,  /* CleanCache 2 */
  };

enum MCOnlineStageType
  { STAGE_BOUND_FIND_UPDATE=0,
    STAGE_UPDATE_BOUND_FIND,
    STAGE_UPDATE_FIND_ALL,
    STAGE_UPDATE_FIND_ONE
  };

typedef struct {
  double beta;
  double epochs;
  long   is_voted;
  KernelDef kernel_def;
  enum RedOptType redopt_type;
  double delta;
  double gamma;

  enum MCOnlineUpdateType  update_type;
  enum MCOnlineStageType   stage_type;
  enum MCOnlineFindType    find_type;
  long spp_pattern_size;
  
} MCOnlineParamDef;


MCSolution mconline(MCDataDef mc_datadef, MCOnlineParamDef mconline_pd, MCClassifier* mc_cls_init);

#endif

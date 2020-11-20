C $Header: /u/gcmpack/MITgcm/pkg/gchem/GCHEM_OPTIONS.h,v 1.9 2011/12/24 01:04:47 jmc Exp $
C $Name:  $

#ifndef GCHEM_OPTIONS_H
#define GCHEM_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

#ifdef ALLOW_GCHEM

CBOP
C    !ROUTINE: GCHEM_OPTIONS.h
C    !INTERFACE:

C    !DESCRIPTION:
c options for biogeochemistry package
CEOP

#define GCHEM_SEPARATE_FORCING

C Preformed tracers typically set to surface level values
C    but here you can set the entire mixed layer.
#undef GCHEM_PREFORMED_MLD
C or just regularly at the surface
#define GCHEM_PREFORMED


#endif /* ALLOW_GCHEM */
#endif /* GCHEM_OPTIONS_H */

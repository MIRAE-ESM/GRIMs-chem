! RTM Check
#if !defined(RRTMGLW) || !defined(RRTMGSW)
#error "Wrong radiation scheme for GRIMs-Chem!"
#error "Use both RRTMGLW and RRTMGSW!"
#endif

! LSM Check
#ifndef NOALSM1
#error "Wrong land scheme for GRIMs-Chem!"
#error "Use Noah land scheme!"
#endif

! CPS Check
#ifndef SAS
#error "Wrong cumulus parameterization scheme for GRIMs-Chem!"
#error "Use SAS scheme!"
#endif

! MPS Check
#if (_nwmass_ < 3)
#error "Wrong microphysics scheme for GRIMs-Chem!"
#error "Use WSM3/WSM5/WSM6/WDM5/WDM6 scheme!"
#endif

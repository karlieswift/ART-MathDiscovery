# Stage05F2C narrow Boltz-2 structure test

Frozen lead: SF_00172800, exact-AA SHA256 73f4e1517b9e8df706aefcf7fa60df9ecda1e0ed210b45eaf84834a0e2e42239, 205 aa, RVT-UG9, E=4.1e-5, qcov=0.5512, model coverage=0.5070, YxDD-negative, replicated in PX404140.1 and PX418225.1.

This package only prepares and predicts one full-3D model under the same primary Boltz-2 inference protocol used in D4-R5: MSA-assisted, seed 20261007, 3 recycles, 200 sampling steps, 1 diffusion sample, `--no_kernels`, GPU. It does not run the calibrated RT-fold decision yet.

Run prepare in `art_rt_evidence`; run Boltz in `art_structure`.

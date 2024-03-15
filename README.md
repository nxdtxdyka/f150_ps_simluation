# f150_ps_simluation
We have simply updated the low-frequency end of the radio continuum simulation in T-RECS (Bonaldi et al., 2018) based on the latest catalog (Best et al., 2023).

### Updated section
We have made simple modifications to T-RECS' random_modules.f90, sampler_continuum.f90, and sampler_modules.f90. Therefore, you only need to perform simple replacements in the src directory of the T-RECS code and then compile it to use.

Specific installation steps can be referred to in the original T-RECS code. The original T-RECS code can be found at the following URL: [[URL]
](https://github.com/abonaldi/TRECS/tree/master)

### Model input files
Model input files, such as the updated luminosity function, can be downloaded from this link: https://www.dropbox.com/scl/fi/i08xe9mziqfjeclipeje9/TRECS_Inputs.7z?rlkey=hhstrpjvzlnuz4psjgzew9c47&dl=0

# The non-is LPJ-GUESS crashes, 11 August 2026

What is written down here is what can be shown from the run directories and the git
metadata. Where something is inference it says so, and where something does not fit it
says that too.

## The runs

| time | job | experiment | signature |
|---|---|---|---|
| 16:39 | 26877117 | 11H0_raupach_null | `Illegal Frac_air -0.142886`, exit 99 |
| 21:08 | 26882104 | 11H0_raupach_null | froze at step 4, killed at 4:17 |
| 21:37 | 26882623 | 11H0_raupach_null | froze at step 4, killed at 3:56 |
| 21:58 | 26883003 | 11H0_codeonly | `Illegal Frac_air -0.245021 / -0.353886`, exit 99 |

These run out of `model_dir: /work/ab0246/a270092/model_codes/awiesm3-develop`, which is
a different tree from the ice sheet work in `awiesm3-develop-is`. Nothing done to the
`-is` tree can reach them.

## What moved, and when

The reflog of `awiesm3-develop/lpj_guess`:

```
07-23 14:43  clone              -> 2e9ce65 (06-15)   June code, no hardening
08-07 16:08  pull: Fast-forward -> 2d01045 (08-07)   brings in cabcaaf and the rest
08-11 13:55  merge              -> d9d1092
```

The 7 August pull is the moment that tree gained every land-fraction hardening commit:
`70c9745`, `561c005`, `7de7005`, `bc36fbd`, `628c2db`, `cabcaaf`, `1cea929`, `2d01045`.

The binaries in that tree do not agree with each other:

```
bin/guess              8c5ab467  07-23 19:00   built from the June code, no hardening
lpj_guess/build/guess  8d2630cc  08-11 21:34   built after the pull, hardened
```

The 11H0 experiment staged `8d2630cc`, in both `bin/lpj_guess/guess` and `work/guess`.
So those runs executed the hardened code against states written by the older code line.
That is the same shape as the ice sheet failure and it is the strongest part of the
case.

## Where the usual description is wrong

The failure has been described as a post-`cabcaaf` check firing. The check is not new.
`Illegal Frac_air` is present in the June code at `2e9ce65`, and `git log -S` traces it
to 10 May and before that to the original 4.1 import in January. The hardening did not
add the abort being hit.

What those commits do is clamp land-cover fractions to `[0,1]` and adjust the tolerance
handling. So the mechanism is more likely that the values arriving at an old check
changed, rather than that a new check appeared. The distinction matters, because it
means reverting to the older binary might not be sufficient, and the states may be
marginal against that check either way.

## What does not fit

The 16:39 failure predates the 21:34 build by five hours. Either an earlier hardened
build existed that cannot now be inspected, since that experiment path was reused by
later attempts and overwritten, or that first failure has a different cause. There is no
way to tell from what survives on disk.

## How this differs from the ice sheet case

In the `-is` setup the pooled state `lpjg_state_3850` fails with a bit-identical value,
`invalid ice volume content -0.016162`, under a binary that predates all the hardening
and under one that includes it. The code is therefore ruled out and the state is at
fault on its own.

Here the failure correlates with the build change and has not been tested against the
older binary. The build is still a live suspect, and the two cases should not be assumed
to have one fix.

## A hazard to be aware of before the next run

`awiesm3-develop` is currently inconsistent: `bin/guess` is the July build and
`lpj_guess/build/guess` is the August one. `bin_dir` is `${setup_dir}/bin`, so which
binary a leg stages depends on whether anyone has copied the build output across. The
build script `comp-lpj_guess-*_script.sh` compiles into `lpj_guess/build/` and does not
install, so a rebuild alone does not reach a run. Check what the leg actually loaded,
with `md5sum <run>/work/guess`, rather than what was built.

## Deliberately excluded

The OOM on 11G leg 4 at 23:49, with `task 2224 oom_kill` and no soil error, and the
cancellation of 11H0_lpjgonly at 00:21 after 80 seconds with no error in any log, are
not included above. Neither carries the soil or land-fraction signature, and there is no
evidence tying them to this. They should be treated separately until something connects
them.

## What would settle it

Run one leg of an affected configuration with the July binary `8c5ab467` and nothing
else changed. If it survives, the build is the cause and the states need respinning
against the hardened code. If it fails the same way, the states are at fault on their
own, as in the ice sheet case, and the build is a red herring.



# sTGC time calibration 

The zero field alignment took place from run `23072003` to `23072023`

## Restore DAQ files

finding the files on hpss:
You can use `hsi` on RCF to get into a terminal for browsing HPSS.

```sh
[rcas6016] ~/> hsi
Username: jdb  UID: 11210  Acct: 11210(11210) Copies: 1 COS: 0 Firewall: off [hsi.8.3.0.p0 Thu Jun 10 07:40:52 EDT 2021]
? 
```

now you are in a terminal for browsing hpss.
the DAQ files are at paths like this:

```sh
ls -la /home/starsink/raw/daq/<year>/<day>/<run_number>/<stream>_<run_number>_<raw|adc>_<evb_id>.daq
```

for instance for run 23072003 the path would be like:

```sh
ls -la /home/starsink/raw/daq/2022/072/23072003/*
```

the evb ids are always the same (unless and evb machine has died - rare). So once you identify the path you can easily get multiple runs, eg.:

```sh
ls -la /home/starsink/raw/daq/2022/072/23072003/st_physics_23072003_raw_1000006.daq

ls -la /home/starsink/raw/daq/2022/072/23072004/st_physics_23072004_raw_1000006.daq

ls -la /home/starsink/raw/daq/2022/072/23072005/st_physics_23072005_raw_1000006.daq
```

Please do not use the daq files from the evb that starts with 00, e.g. st_physics_23072003_raw_0000003.daq, since those are usually files with few events.

finally, to restore them use the `hpss_user.pl` command (on RCF) :

```sh
[rcas6016] ~/> hpss_user.pl <full_path_on_hpss> <full_path_on_RCF>
```

For instance:

```sh
hpss_user.pl /home/starsink/raw/daq/2022/072/23072003/st_physics_23072003_raw_1000006.daq /star/data03/pwg/jdb/FWD/daq/st_physics_23072003_raw_1000006.daq
```

then you can check progress with:

```sh
hpss_user.pl -w


Date > May 10 2022 and Expiration > May 18 2022

                   Status       User   Count  Sanity
 ---------------------------------------------------------------------------
                      New         li       1  New and never submitted
                             starlib    7986  New and never submitted
                             starlib    3650  Requeued
                             starlib    6183  Requeued after 1 failure
                            starreco    2267  New and never submitted
 ---------------------------------------------------------------------------
    Submitted/In Progress    starlib     251  Waiting for HPSS response (new mode)
                            starreco      50  Waiting for HPSS response (new mode)
 ---------------------------------------------------------------------------
                 Restored        jdb       4  Success on first try (new mode)
                               nihar       4  Success on first try (new mode)
                             starlib  149426  Success on first try (new mode)
                             starlib    8261  Previously failed 1 time
                            starreco   52970  Success on first try (new mode)
                            starreco    4474  Previously failed 1 time
                             tinglin    1116  Success on first try (new mode)
                              tlusty     132  Success on first try (new mode)
                              tlusty      41  Previously failed 1 time
                               vipul       8  Success on first try (new mode)
                               vipul       7  Previously failed 1 time
                                 yxu     426  Previously failed 1 time
 ---------------------------------------------------------------------------
                  Failure    saskiam      10  Previously failed 1 time
                          sethcarl20      10  Non-recoverable error occured (new mode)
                             starlib    6924  Non-recoverable error occured (new mode)
                             starlib    1183  Previously failed 1 time
                            starreco       1  Non-recoverable error occured (new mode)
                            starreco      14  Previously failed 1 time
                              tlusty       3  Previously failed 1 time
 ---------------------------------------------------------------------------
       Marked for requeue    starlib    1756  Requeued known recoverable (new mode)
                             starlib     313  Previously failed 1 time
 ---------------------------------------------------------------------------
      Stalled or vanished      nihar       1  Previously submitted 1 time
                             starlib   16175  Previously submitted 2 times
                            starreco     644  Previously submitted 1 time
                              tlusty       8  Previously submitted 1 time
                                 yxu      89  Previously submitted 2 times
 ---------------------------------------------------------------------------
                    Bogus      nihar       5  Previously submitted 2 times
                          sethcarl20      20  Previously submitted 1 time
                               vipul       1  Previously submitted 1 time
-
                      New %tage=07.59%
    Submitted/In Progress %tage=00.11%
                 Restored %tage=82.01%
                  Failure %tage=03.08%
       Marked for requeue %tage=00.78%
      Stalled or vanished %tage=06.39%
                    Bogus %tage=00.00%
```


Scratch is available on `/star/data05/scratch/`
just make your own sub directory there, e.g.:

```sh
mkdir /star/data05/scratch/stgc_zf_calib
```

then use that directory for restoring the DAQ files.

For the time calibration we only need 2 files from every run between `23072003` to `23072023`.
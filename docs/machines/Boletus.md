# Boletus

These are detailed information about each available machine. **DO NOT** put passwords, keypair locations or other sensitive information on these pages - they are public!

## Basic info

Name on network | IP  | Name | Total RAM(GB) | Total cores | Notes
--------------- | --- | ---- | ------------- | ----------- | -----
boletus.ad.kew.org | x.x.1.102 | Boletus | 256 | 24 | na

## Known Sudoers

Joe, Alex, Mike, Pepijn

## Specific notes for this machine

(Put any notes specific to this machine here.)

* Update 2016-aug-19*: installed poretools 0.5.1 from github (Joe)
* Update 2016-dec-14: installed OpenCL using apt-get, CUDA 8.0 following instructions on https://developer.nvidia.com/cuda-downloads (Pepijn)
* Update 2016-dec-15: installed parallel version of MrBayes 3.2.6 -> type `mpirun -np NUMBER_CPUs mbpar` (Pepijn)
* Update 2017-jan-04: installed tuptime 3.3.0 to monitor shutdowns/reboots -> type `sudo tuptime` (Pepijn)
* Update 2017-mar-21: updated boost to v 1.63.0, updated mothur to v 1.39.5, updated uchime to v 4.2, updated FastQC to v 0.11.5 (Pepijn)
* Update 2017-mar-24: updated ABySS to v 2.0.2, available from /usr/local/bin -> version in /usr/bin is still v 1.5.2 (Pepijn)
* Update 2017-jun-29: deployed Trinity from github in /usr/local/bioinf/git-repositories/, available as symlink trinity from /usr/local/bin (Pepijn)
* Update 2017-jul-31: upgraded ant to v 1.10.1 and installed jModelTest2 from github in /usr/local/bioinf/git-repositories/, available as symlink jModelTest2 (Pepijn)
* Update 2017-nov-13: installed nanopolish from github in /usr/local/bioinf/git-repositories/, available as symlink nanopolish (Pepijn)
* Update 2017-nov-14: installed HTSlib v1.6.24 and samtools v1.6.2 from github in /usr/local/bioinf/git-repositories/ (make install into /usr/local/bin/), to work with nanopolish. Old version samtools available as samtools_0.1.19 (Pepijn)

**Student-preference machine**: This machine is primarily indended for MSc students' projects. In all cases they are assumed to have priority in job scheduling - please check before running long or CPU- or memory-intensive analyses.


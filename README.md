# FlexAID
Hi NRG’ers,

The latest source code of FlexAID∆S is available on the NRGlab Github under the Entropy branch (I kept the master branch as FlexAID published version by Francis). You can clone the whole FlexAID∆S repertory using the following command

    > git clone https://github.com/NRGlab/FlexAID.git

 and switching to branch Entropy using 

    > git pull Entropy && git checkout Entropy

Alternatively you can download it as a zip @ : https://github.com/NRGlab/FlexAID/archive/Entropy.zip

The problem was with the compilation of the original FlexAID and the makefile didn’t compile FlexAID without Entropy files. I pushed commits to the master (FlexAID) and the Entropy branch (FlexAID∆S) Friday and it seems to compile fine everywhere but ComputeCanada using the Makefile I pushed Friday. Boost, our only dep, has changed since I last checked and our Makefile is outdated.

Use the command `module load boost/1.68.0.lua` on Compute Canada servers exemple Beluga. Tested it last week, should be fine.

Compute Canada clusters uses the command the command 

`module load boost/1.68.0.lua` and it installs in standard locations to ease compilation.

You can compile FlexAID by changing directory into the BIN subdirectory and using the following command:

`make -f Makefile.Linux64 clean && make -f Makefile.Linux64`

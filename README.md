
**This repository has been included in the CMSHI repository.**

**https://github.com/CmsHI/RawToMP7Text.**

**Please refer to this repository for the current version of this code**

# Instructions ConvertedoRSCtoMP7


This code allows you to read an oRSC output txt file and convert it into the MP7 payload format.

The basic class is DefineLinks.cc.


In InputsoRSC you can find some examples of oRSC input files for different MC pp events 
produced with the UCT RCT emulator. 
In the OutputMP7 folder you will find the converted files in the MP7 format.

In the OldCode folder, I stored a previous version of the DefineLinks class.

Instructions to run the code:

 * First clone the repository with: 
          git clone https://github.com/ginnocen/ConvertedoRSCtoMP7
 * Compile the code
          cd test
          make
 * You can run the executable file:
          ./testGlobal

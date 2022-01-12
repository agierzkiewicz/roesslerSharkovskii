# A. Gierzkiewicz, P. Zgliczyński, "From the Sharkovskii theorem to periodic orbits for the Rössler system"

## An installation instruction of the programs realizing the numerical part of the proofs of contracting grids' and periodic orbits' existence on a Poincaré section.

### Requirements:
The program is written in C++ and has been tested under Linux Mint 18.1 with gcc (Ubuntu 5.4.0-6ubuntu1~16.04.12) compiler. The program uses the [CAPD](http://capd.ii.uj.edu.pl/) library ver. 5.0.6 (see also [sourceforge](https://sourceforge.net/projects/capd/files/) download zone). The CAPD library is also available as debian deb package. There are also prebuilt versions for Debian, Ubuntu and OSX.

### Compilation instruction:
In principle it is possible to compile the CAPD library under MS Windows but we strongly recommend to compile and run the programs on linux-like systems. The following commands should be executed from the terminal (or msys environment under MS Windows)

- unzip the archive: ```unzip Roessler_Sharkovskii.zip```
- enter directory of the program: ```cd Roessler_Sharkovskii```
- build the programs: ```make CAPDBINDIR=relative_or_absolute_path_to_capd_bin_directory```

    for example: ```make CAPDBINDIR=../../capd-build/bin/```
    
    Warning: Do not forget the last slash character.
- If the CAPD library is installed in standard directories that are on system path (like /usr/bin/) then you can compile the program just by invoking ```make```
    
    The command generates executable files in the current directory.

### Program options:
There is one program with computer-assisted proofs in the package: ```Roessler_Sharkovski.cpp```. 
It can be run from the terminal by ```./Roessler_Sharkovskii```.

The main menu asks to choose one of six procedures:

    1. Proof of the contracting grids's existence for a = 5.25
    2. Proof of the contracting grids's existence for a = 4.7
    3. Proof of the 6-periodic orbit's existence by INM for a = 4.381
    4. Proof of the contracting grids's existence for a = 4.381
    5. Proof of the 6-periodic orbit's existence by INM for a = 5.42
    6. Proof of the contracting grids's existence for a = 5.42 

### Output:

The program executes within 6 minutes on a laptop type computer with Intel Core i7 2.7GHz x 2 processor (Procedure 2 takes 1:50 min and Procedure 4 takes 3:30 min).

The program does not create output files. It prints to the screen what conditions are being checked and the result of this check (true or false).

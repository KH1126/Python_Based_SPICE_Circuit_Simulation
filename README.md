# Python_Based_SPICE_Circuit_Simulation

All source code are in src folder

I.AC analysis :
1. run test_func.py
2. enter circuit file's name
3. enter command as following template
  ac <variation type> <number of points> <start frequency> <stop frequency>
  Eg. ac dec/oct/lin 100 1 100000
notice : frequency must > 0

II. DC analysis
1. run test_func.py
2. enter circuit file's name
3. enter command as following template
  dc <source> <start value> <stop value> <increment>
  Eg. dc Vin 0 2 0.01
notice : If end voltage<start voltage increment must < 0
  Eg. dc Vin 2 0 -0.01

III. Transient analysis
1. run transient_0327.py
2. enter circuit file's name
3. enter command as following template  tran <step size> <stop time>
   Eg. tran 0.1 1
   notice : The unit above is sec.
   

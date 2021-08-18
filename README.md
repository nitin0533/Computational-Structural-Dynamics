# Computational-Structural-Dynamics
## Resume points
**Component Mode Synthesis to estimate Natural frequencies of Structural Assembly using MATLAB**
1. Performed Model Order Reduction of assembly by 18 times using Craig-Bampton Method & found First 10 Frequencies
2. Improved CBM results with Model Order reduction of assembly by 22 times using Characteristic Constraint Modes
3. Acheived 98.34% accuracy in first 10 Natural frequencies using CBM, and improved it to 99.97% using CC modes

**Component Mode Synthesis to estimate Natural frequencies of a Structural Assembly using MATLAB**
1. Performed Model Order Reduction by 18 times using Craig-Bampton Method and 22 times using Characteristic Constraint modes
2. Acheived 98.34% accuracy in first 10 Natural frequencies using CBM, and improved it to 99.97% using CC modes
 ![mass_proof](https://user-images.githubusercontent.com/71177034/129353591-ec4824b0-7ef8-49f2-bb25-4bc0600b434d.png)
  ![stiffness_proof](https://user-images.githubusercontent.com/71177034/129355195-ae507b19-232f-4327-b865-215a44bd3ad5.png)
![stiffness_proof_2](https://user-images.githubusercontent.com/71177034/129355214-9a438fca-b4cc-4a72-a559-848645a2f869.png)



Screenshots of the reports and MATLAB workspace when the code **19310R001_code.m** is run is attached. 
 
 M and K --> mass and stiffness matrices of the assembled structural system (order=3062x3062)
 
 M_cb and K_cb --> mass and stiffness matrices of the structural system by Craig-Bampton method (order=166x166)
 
 M_cc and K_cc --. mass and stiffness matrices by Charateristic constraint (CC) modes (order=140x140)
 
 By Craig-Bampton method the order of matrices is reduced by (3062/166=18.44 rounded to 18) times
 
 By Characteristic constraint modes the order of matrices is reduced by (3062/140=21.87 rounded to 22) times
 

 ![proof](https://user-images.githubusercontent.com/71177034/129355356-92a7016e-ab2b-4e96-a4d7-a19e5daf646b.png)


The results of first ten natural frequencies are tabulated below (snip from the report)
![image](https://user-images.githubusercontent.com/71177034/129355579-0da2d88a-8161-45bd-b4ff-1a9e07fb6e80.png)


Minimum accuracy between Craig-Bampton frequencies and actual frequencies = (1-(32.9365-32.3998)/32.3998)*100 =98.34%

Minimum accuracy between CC modes frequencies and actual frequencies = (1-(32.4003-32.3998)/32.3998)*100 = 99.97%

**How to run the code**
 Make sure all the files are in same directory. Run the code.

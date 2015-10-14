RPG README
May 5 2010
---

Guidelines for humans model generation

1) Read the Humans documentation 'UserDocumentation.pdf' 
chapter 2 (p21-32) Computation of Lagrangian Model
The directory  MapleCodeGeneration/ is a copy of the Humans
one in humans/LagrangianModel/PA10/MapleCodeGeneration/

2) Files in MapleCodeGeneration/

input maple files:
 KinematicData.maple - define the kinematic of the robot = Geometric modelisation of a jointed system (see chapter2 for the parametrization)

 AdditionnalData.maple - characteristic points, called tags, were attached on segment of the model; There are used for visualization or define contacts.

 DynamicData.maple - define the kinematic of the robot = In this file, the inertials parameters (the gravity vector, the mass of the segment, the position of the segments centers of mass and the inertia matrix of these segments relative to the centers of their attached frames) are defined. 

maple generators files:

output C files:
Contact.c         
ContactJacobian.c    
NLEffects.c
ContactHessian.c  
HAnimOrientations.c 
Inertia.c        
Tags.c


3) How to generate the code

* erase or rename  output C files.

* xmaple GenerateModel.mw
  in the gui, relaunch 
  restart:
  read `KinematicData.maple`;
  read `DynamicData.maple`;
  read `AdditionnalData.maple`;
  read `ModelGeneration.maple`;
  GenerateModel();



#!/usr/bin/env python


import numpy
import array
my_PI=3.14159265

#BODIES DESCRIPTION
NBBODIES=8

BATIE=-1
MANETTE=0
BIELLETTE=1
CROCHET=2
PLATINE=3
EQUIPAGE=4
BARRE=5
PROBE=6
AIGUILLE=7


body=numpy.array(['manette','biellette','crochet','platine','equipage','barre','probe', 'aiguille'])
##### Initial body positions in  ON CONDITION ########################################################
initPos=numpy.array([(-16.1,10,0.5,0,0,1,-19.0*my_PI/180.0),
                     (-11.659106,8.461021,10.3,   0,0,1, -5.47692     *my_PI/180.0),
                     (7.766977,7.637695,7.1,     0,0,1, 14.6953    *my_PI/180.0),
                     (0,0,9.7,           0,0,1, -2        *my_PI/180.0),
                     (0,0,0,             0,0,1, -71.25       *my_PI/180.0),
                     (7,-0.25,7.1,       0,0,1, 36        *my_PI/180.0),
		     (-5.9,-8.9,4.65,           0,1,0, -180*my_PI/180.0),#-2.429
		     (0,0,5.1,           1,0,0, 90*my_PI/180.0),
])

initVel=numpy.array([(0,0,0,0,0,0),
		     (0,0,0,0,0,0),
		     (0,0,0,0,0,0),
		     (0,0,0,0,0,0),
		     (0,0,0,0,0,0),
		     (0,0,0,0,0,0),
                     (0.0,0,0,0,0,0), ## ON-OFF => -2 ## OFF-ON =>2
		     (0,0,0,0,0,0),

])
###### CENTRE OF MASS of bodies  ########################################################
initCenterMass=numpy.array([(1.123281e+00, 2.939292e+00, 4.001829e+00),
                            (7.562450e+00, -1.224523e+00, -1.585747e-01),
                            (-1.966777e+00, 7.198357e-01, 1.282473e+00),
                            (3.319088e+00, 4.365461e+00, -4.706395e+00),
                            (7.7911593e+00, -1.7807128e+00,  4.0343125e+00),
                            (-1.9723481e+00, -2.3250491e+00, -4.9871647e+00),
			    (2.844,0,0),
			    (0,0,5),
			    #(0.0000000e+00,  -8.0000000e+00,   0.0000000e+00)
			    #(0,0,5)
])

##### MASS of bodies ###############################################################
m=numpy.array([1.30,
               0.50,
               0.09,
               1.35,
               2.50,
               0.77,
	       0.8,
	       0.1,
])
###### INERTIA MATRIX of bodies in (Ixx,Ixy,Ixz),(Iyx,Iyy,Iyz),(Izx,Izy,Izz) #################################################
inertialMatrix=numpy.array([((3.891182e+01, -6.015089e+00, -2.8917964e-01),(-6.0150893e+00,  2.7959058e+01, -1.9219978e+00),( -2.8917964e-01, -1.9219978e+00,  4.1931367e+01)),
                            ((2.4550110e+00,  3.1330348e-01, -1.4650670e+00),(3.1330348e-01,  2.1629420e+01,  6.8441857e-01),(-1.4650670e+00,  6.8441857e-01,  2.0917600e+01)),
                            ((3.0507502e-01,  8.1781430e-02, -5.4128496e-02),(8.1781430e-02,  4.4192990e-01,  1.9009435e-02),(-5.4128496e-02,  1.9009435e-02,  6.2063100e-01)),
			    ((2.9772412e+01,  3.5465101e+00,  4.7496642e-01),(3.5465101e+00,  2.5095820e+01,  6.0505649e-01),(4.7496642e-01,  6.0505649e-01,  3.7698441e+01)),
                            ((2.3767882e+01,  8.5414908e-01, -1.1814918e+01),(8.5414908e-01,  1.7702250e+02, -1.4361701e-01),(-1.1814918e+01, -1.4361701e-01,  1.8255743e+02)),
                            ((18.429238, -1.9231011,  5.9362999),(-1.9231011,  31.822293, -5.3010733e-01),(5.9362999, -5.3010733e-01,  2.7636361e+1)),
			    ((1.78e+01,0,0),(0,1.61e+02,0),(0,0,1.61e+02)),
			    ((1,0,0),(0,1,0),(0,0,1)),
			    #((8.3896687e+00,  0.0000000e+00, -8.7749247e-6),(0.0000000e+00,  1.9192988e-01,  0.0000000e+00),(-8.7749247e-6,  0.0000000e+00,  8.3896744e+00))
])
###### STEP Files of bodies ################################################################
afile=['./CAD/C60/manette.stp',
                   './CAD/C60/biellette.stp',
                   './CAD/C60/crochet.stp',
                   './CAD/C60/platine.stp',
                   './CAD/C60/portcontact/porte_contact.stp',
                   './CAD/C60/newFC/Tbarassembly.stp',
		   './CAD/C60/C8/probe.stp',
		   './CAD/C60/J1/Aiguille.stp',
		   #'./CAD/C60/J1/72569701.stp'
#		   './CAD/C60/J1/couvercle.stp' # porte_contact.stp=> (0,0,0,0,0,1, -69*my_PI/180.0),
]
###### FORCES and MOMENTS acting on the bodies##################################################
plugin='plugins.so'

fctfext=numpy.array(['',
                  '',
                  '',
                  '',
                  '',
                  '',
		  '',
		  ''#
])

fctmext=numpy.array(['',
                  '',
                  '',
                  '',
                  '',
                  '',
                  '',
		  ''
])

fctfint=numpy.array(['externalManetteForces',
                     '',
                     '',
                     'externalPlatineForces',
                     'externalEquipageForces',
                     'externalBarreForces',
                     '',
                     ''#
])

fctmint=numpy.array(['externalManetteMomentum',
                     '',
                     '',
                     'externalPlatineMomentum',
                     '',
                     'externalBarreMomentum',
		     '',
                     ''
                 ])

fctfintjacq=numpy.array(['FiniteDifference',
                         '',
                         '',
                         'FiniteDifference',
                         'FiniteDifference',
                         'FiniteDifference',
                         '',
                         ''
                         ])
fctmintjacq=numpy.array(['FiniteDifference',
                         '',
                         '',
                         'FiniteDifference',
                         '',
                         'FiniteDifference',
                         '',
                         ''])

## REQUIRED Boundary condition
# we impose a boundary condition for the foirst body given by the plugin
boundaryCondition=numpy.array(['','','','','','','prescribedvelocityB1',''])
# we prescribe the 4th component of the velocity of the first body
# i.e we prescribe a angular velocity around the y-axis.
boundaryConditionIndex=numpy.array([numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([]),numpy.array([0]), numpy.array([])]) # we impose the velocity

###### VISIBILITY of CAD geometry###########################################################
bodyDraw=array.array('I',[1,
                          0,
                          1,
                          1,
                          1,
                          1,
                          1,
			  1
])
                  
######JOINT DESCRIPTIONS#################################################################
#JOINT DESCRIPTION
NBJOINTS=0
jointName=numpy.array([#'Manette0',
                       #'ManetteBiellette',
                       #'BielletteCrochet',
                       #'Platine0',
                       'Equipage0',
                      #'CrochetPlatine',
                       #'PlatineBarre'
])

jointType=array.array('I',[#mbtb.PIVOT_0,
                       #mbtb.PIVOT_1,
                       #mbtb.PIVOT_1,
                       #mbtb.PIVOT_0,
                       mbtb.PIVOT_0,
                       #mbtb.PIVOT_1,
                       #mbtb.PIVOT_1
])

jointPos=numpy.array([#(0,0,1,0,0,0),
                      #(0,0,1,  4.7,0,0),
                      #(0,0,1,   16.6,0,0),
                      #(0,0,1,   0,0,0),
                      (0,0,1,   0,0,0),
                      #(0,0,1,   0,0,0),
                      #(0,0,1,   7,0,0)
])
jointBody1=array.array('I',[#MANETTE,
                            #MANETTE,
                            #BIELLETTE,
                           # PLATINE,
                            EQUIPAGE,
                            #CROCHET,
                            #PLATINE
])

jointBody2=array.array('I',[#0,
                            #BIELLETTE,
                            #CROCHET,
                            #0,
                            0,
                            #PLATINE,
                            #BARRE
])

######CONTACTS DESCRIPTIONS#################################################################
#CONTACTS DESCRIPTION
NBCONTACTS=39 #40/23/38

contactName=numpy.array([
	'Handle_case_J1_bot',
	'Handle_cover_J1_top',
	'Case_handle_C1_ON_1',
	'Case_handle_C1_ON_2',
	'Cover_handle_C1_ON_3',
	'Cover_handle_C1_ON_4',
	'Case_needle_S3_bot',
	'Cover_needle_S3_top',
	'Needle_plate_J2_bot',
	'Needle_plate_J2_top',
	'Plate_hook_J3_bot',
	'Plate_hook_J3_top',
	'Plate_Tbar_J4_bot',
	'Plate_Tbar_J4_top',
	'Hook_Tbar_C4_1',
	'Hook_Tbar_C4_2',
	'Hook_Tbar_C4_3',
	'Hook_Tbar_C4_4',
	'Rod_handle_J5_bot',
	'Rod_handle_J5_top',
	'Rod_hook_J6_bot',
	'Rod_hook_J6_top',
	'Mcontact_needle_J7_bot',
	'Mcontact_needle_J7_top',
	'Mcontact_Fcontact_bot_1_C7',
	'Plate_case_C3_ON_1',
	'Plate_case_C3_ON_2',
	'Plate_case_C3_ON_3',
	'Plate_case_C3_ON_4',
	'Tbar_Hook_C5_1',
	'Tbar_Hook_C5_2',
	'Mcontact_plate_C6_side_1',
	'Mcontact_plate_C6_side_2',
	'Mcontact_plate_C6_side_3',
	'Mcontact_plate_C6_side_4',
	#'Probe_FixedContact',
	#'Tbar_FC_X',
	#'Tbar_FC_Y',
	#'Tbar_Lamage_FC_Y',
	#'Tbar_Magnetic_FC_X',
	#'Mcontact_Magnetic_FC_X',
	#'Tbar_Lamage_FC_axis_X',
	#'Tbar_Lamage_FC_axis_Y',
	#'Mcontact_Case_stopper1', #S71#
	#'Mcontact_Case_stopper2', #S72#
	##'Rod_Cover_top1', #S6
	##'Rod_Cover_top2', #S6
	#'Plate_Cover_top',#S11
	#'Tbar_Cover_top',#S12
	#'Handle_Case_bottom_surface', #S1
	#'Hook_Plate_bottom_surface',#S41
	#'Hook_Plate_bottom_surface',#S42
	#'Tbar_Plate_bottom_surface',#S51
	#'Tbar_Plate_bottom_surface',#S52
	#'Plate_Mcontact_bottom_surface',#S81
	#'Plate_Mcontact_bottom_surface',#S82
	#'Hook_Cover_top_surface',#Sx
	'TRIP_1',
	'TRIP_2',
	'Tbar_Magnetic_FC_X',
	'plunger_bearing'



])
######CONTACT BODY1 STEP FILES DESCRIPTIONS#################################################################
afileContact1=[
	'./CAD/REVISED_CAD2/J1/Handle_Case_J1_bottom_bearing2.stp',
	'./CAD/REVISED_CAD2/J1/Handle_Case_J1_top_bearing2.stp',
	'./CAD/C60/C1/Handle_case_C1_ON_plane.stp',
	'./CAD/C60/C1/Handle_case_C1_ON_plane.stp',
	'./CAD/manette_handle_C2_1.stp',
	'./CAD/manette_handle_C2_2.stp',
	'./CAD/C60/S3/Needle_journal_S3_surface.stp',
	'./CAD/C60/S3/Needle_journal_S3_surface.stp',
	'./CAD/REVISED_CAD2/J2/Plate_Needle_J2_bottom_braring2.stp',
	'./CAD/REVISED_CAD2/J2/Plate_Needle_J2_top_braring2.stp',
	'./CAD/REVISED_CAD2/J3/Hook_J3_bottom_ring2.stp',
	'./CAD/REVISED_CAD2/J3/Hook_J3_top_ring2.stp',
	'./CAD/REVISED_CAD2/J4/Tbar_J4_bottom_ring2.stp',
	'./CAD/REVISED_CAD2/J4/Tbar_J4_top_ring2.stp',
	'./CAD/REVISED_CAD3/C4/Hook_half_ring_C4_bottom.stp',#/Hook_latest/Hook_Tbar_bottom_ring_1.stp',#
	'./CAD/REVISED_CAD3/C4/Hook_half_ring_C4_top.stp',#/C60/Hook_latest/Hook_Tbar_bottom_ring_2.stp',#
	'./CAD/REVISED_CAD1/C4/Hook_Tbar_C4_surface.stp',
	'./CAD/REVISED_CAD1/C4/Hook_Tbar_C4_surface.stp',
	'./CAD/REVISED_CAD2/J5/Handle_J5_bottom_bearing2.stp',
	'./CAD/REVISED_CAD2/J5/Handle_J5_top_bearing2.stp',
	'./CAD/REVISED_CAD2/J6/Hook_Rod_J6_bottom_bearing2.stp',
	'./CAD/REVISED_CAD2/J6/Hook_Rod_J6_top_bearing2.stp',
	'./CAD/REVISED_CAD2/J7/Mcontact_Needle_J7_bottom_bearing2.stp',
	'./CAD/REVISED_CAD2/J7/Mcontact_Needle_J7_top_bearing2.stp',
	'./CAD/C60/portcontact/porte_contact_asm_contact_surface.stp',
	'./CAD/C60/C3/Plate_case_stopper_plane_C3.stp',
	'./CAD/C60/C3/Plate_case_stopper_plane_C3.stp',
	'./CAD/C60/C3/Plate_case_stopper_inside_ring_C3.stp',
	'./CAD/C60/C3/Plate_case_stopper_outside1_ring_C3.stp',
	'./CAD/REVISED_CAD4/C5/Hook_Tbar_surface_C5.stp',
	'./CAD/REVISED_CAD4/C5/Hook_Tbar_surface_C5.stp',
	'./CAD/C60/portcontact/porte_contact_asm_stopper_curve1.stp', #C6#
	'./CAD/C60/portcontact/porte_contact_asm_stopper_curve2.stp', #C6#
	'./CAD/C60/portcontact/porte_contact_asm_stopper.stp', #C6#
	'./CAD/C60/portcontact/porte_contact_asm_stopper.stp', #C6#
	##'./CAD/72333401.stp',
	#'./CAD/C60/FC/Tbar_FC_measurement.stp',
	#'./CAD/C60/FC/Tbar_FC_measurement.stp',
	#'./CAD/C60/newFC/Tbar_lamage_FC_Y.stp',
	#'./CAD/C60/newFC/Tbar_magnetic_FC_X.stp',
	#'./CAD/C60/newFC/porte_contact_FC_X.stp',
	#'./CAD/C60/newFC/tbarassembly_asm_axis_lamage.stp',
	#'./CAD/C60/newFC/tbarassembly_asm_axis_lamage.stp',
	#'./CAD/C60/S7/Movingcontact_case_surfacecontact_ring1_S7.stp',#S71
	#'./CAD/C60/S7/Movingcontact_case_surfacecontact_ring2_S7.stp',#S72
	##'./CAD/C60/S6/Rod_Cover_contact_surface1_S6.stp',#S6
	##'./CAD/C60/S6/Rod_Cover_contact_surface2_S6.stp',#S6
	#'./CAD/C60/S11/Plate_top_surface_cover_side_contact.stp',#S11
	#'./CAD/C60/S12/Tbar_top_surface_cover_side_contact.stp',#S12
	#'./CAD/C60/S1/Handle_bottom_case_S1.stp',#S1
	#'./CAD/C60/S4/Hook_plate_bottomsurface_contact1_S4.stp',#S41
	#'./CAD/C60/S4/Hook_plate_bottomsurface_contact2_S4.stp',#S42
	#'./CAD/C60/S5/Tbar_plate_bottomcontact1_S5.stp',#S51
	#'./CAD/C60/S5/Tbar_plate_bottomcontact2_S5.stp',#S52
	#'./CAD/C60/S8/Plate_movingcontact_surface_contact1_S8.stp',#S81
	#'./CAD/C60/S8/Plate_movingcontact_surface_contact1_S8.stp',#S82
	#'./CAD/C60/Sx/Hook_top_surface.stp',#Sx
	'./CAD/C60/C8/barre_magnetic_3.stp',
	'./CAD/C60/C8/barre_magnetic_2.stp',
	'./CAD/C60/newFC/Tbar_magnetic_FC_X.stp',
	'./CAD/C60/Magnetic_plunger_outer_casing/magnetic_plunder_surface.stp',
]
######CONTACT BODY2 STEP FILES DESCRIPTIONS#################################################################
afileContact2=[
	'./CAD/C60/J1/Case_handle_J1_journal_surface.stp',
	'./CAD/C60/J1/Case_as_cover_handle_journal_surface_J1.stp',
	'./CAD/Case_C2_ring_3.stp',
	'./CAD/Case_C2_ring_4.stp',
	'./CAD/REVISED_CAD/C2/Cover_as_Case_C2_stopper.stp',
	'./CAD/REVISED_CAD/C2/Cover_as_Case_C2_stopper.stp',
	'./CAD/C60/S3/Case_needle_S3_bearing.stp',
	'./CAD/C60/S3/boitier_cover_top_aiguille_surface.stp',
	'./CAD/C60/J2/Needle_journal_S3_surface.stp',
	'./CAD/C60/J2/Needle_journal_S3_surface.stp',
	##'./CAD/Case/case_needle_aasembly_needle_surface.stp',
	##'./CAD/Case/case_needle_aasembly_needle_surface.stp',
	'./CAD/REVISED_CAD/J3/Plate_Hook_J3_journal_surface.stp',
	'./CAD/REVISED_CAD/J3/Plate_Hook_J3_journal_surface.stp',
	'./CAD/REVISED_CAD/J4/Plate_Tbar_J4_journal_surface.stp',
	'./CAD/REVISED_CAD/J4/Plate_Tbar_J4_journal_surface.stp',
	'./CAD/C60/C4/Tbar_hook_wedge_surface_C4.stp',
	'./CAD/C60/C4/Tbar_hook_wedge_surface_C4.stp',
	'./CAD/REVISED_CAD3/C4/Tbar_half_ring_C4_bottom.stp',
	'./CAD/REVISED_CAD3/C4/Tbar_half_ring_C4_top.stp',
	'./CAD/REVISED_CAD/J5/Rod_Handle_J5_Journal_surface.stp',
	'./CAD/REVISED_CAD/J5/Rod_Handle_J5_Journal_surface.stp',
	'./CAD/REVISED_CAD/J6/Rod_Hook_J6_journal_surface.stp',
	'./CAD/REVISED_CAD/J6/Rod_Hook_J6_journal_surface.stp',
	'./CAD/C60/J7/Needle_journal_S3_surface.stp', 
	'./CAD/C60/J7/Needle_journal_S3_surface.stp',
	##'./CAD/Case/case_needle_aasembly_needle_surface.stp',
	##'./CAD/Case/case_needle_aasembly_needle_surface.stp',
	'./CAD/C60/C7/Case_as_fixedcontact_movingcontact_stopper_C7.stp',
	'./CAD/C60/C3/Case_plate_inside_topring_C3.stp',
	'./CAD/C60/C3/Case_plate_edge_topring_C3.stp',
	'./CAD/C60/C3/Case_plate_stopper_surface_C3.stp',
	'./CAD/C60/C3/Case_plate_stopper_surface_C3.stp',
	'./CAD/REVISED_CAD4/C5/Tbar_Hook_half_circular_ring_C5_bottom_side3.stp',
	'./CAD/REVISED_CAD4/C5/Tbar_Hook_half_circular_ring_C5_top_side3.stp',
	'./CAD/C60/C6/Plate_movingcontact_sidesurface_stopper_C6.stp',
	'./CAD/C60/C6/Plate_movingcontact_sidesurface_stopper_C6.stp',
	'./CAD/C60/C6/Plate_movingcontact_sidestop_bottom_ring_C6.stp',
	'./CAD/C60/C6/Plate_movingcontact_sidestop_top_ring_C6.stp',
	##'./CAD/C60/C7/Case_as_fixedcontact_movingcontact_stopper_C7.stp',
	#'./CAD/C60/newFC/Case_Tbar_FC_X.stp',
	#'./CAD/C60/newFC/Case_Tbar_FC_Y.stp',
	#'./CAD/C60/newFC/Case_Tbar_FC_Y.stp',
	#'./CAD/C60/newFC/Case_Magnetic_Mcontact_FC_X.stp',
	#'./CAD/C60/newFC/Case_Magnetic_Mcontact_FC_X.stp',
	#'./CAD/C60/newFC/Case_Tbar_FC_X.stp',
	#'./CAD/C60/newFC/Case_Tbar_FC_Y.stp',
	#'./CAD/C60/S7/Case_movingcontact_bottom_surface_plane_S7.stp',#S71
	#'./CAD/C60/S7/Case_movingcontact_bottom_surface_plane_S7.stp',#S72
	##'./CAD/C60/S6/Case_as_Cover_Rod_Top_surface_stopper.stp',#S6
	##'./CAD/C60/S6/Case_as_Cover_Rod_Top_surface_stopper.stp',#S6
	#'./CAD/C60/S11/Case_as_Cover_Plate_Top_surface_stopper.stp',#S11
	#'./CAD/C60/S12/Case_as_Cover_Tbar_Top_surface_stopper.stp',#S12
	#'./CAD/C60/S1/Case_handle_bottom_stopper_surface_S1.stp',#S1
	#'./CAD/C60/S4/Plate_hook_bottom_surface_S4.stp',#S4
	#'./CAD/C60/S4/Plate_hook_bottom_surface_S4.stp',#S4
	#'./CAD/C60/S5/Plate_Tbar_bottom_surface_S5.stp',#S5
	#'./CAD/C60/S5/Plate_Tbar_bottom_surface_S5.stp',#S5
	#'./CAD/C60/S8/Movingcontact_plate_surfacecontact_S8.stp',#S8
	#'./CAD/C60/S8/Movingcontact_plate_surfacecontact_S8.stp',#S8
	#'./CAD/C60/Sx/Case_Hook_top_surface.stp',#Sx
	'./CAD/C60/C8/probe_surface.stp',
	'./CAD/C60/C8/probe_surface.stp',
	'./CAD/C60/newFC/Case_Magnetic_Mcontact_FC_X.stp',
	'./CAD/C60/Magnetic_plunger_outer_casing/case_plunger_cover.stp',

]
######NAME OF CONTACT BODY1  DESCRIPTIONS#################################################################
contactBody1=array.array('I',[
	MANETTE, #J1
	MANETTE, #J1
	MANETTE, #C1
	MANETTE, #C1
	MANETTE, #C1
	MANETTE, #C1
	AIGUILLE, #S3
	AIGUILLE, #S3
	PLATINE, #J2
	PLATINE, #J2
	CROCHET, #J3
	CROCHET, #J3
	BARRE, #J4
	BARRE, #J4
	CROCHET, #C4
	CROCHET, #C4
	CROCHET, #C4
	CROCHET, #C4
	MANETTE,#J5
	MANETTE,#J5
	CROCHET,#J6
	CROCHET,#J6
	EQUIPAGE, #J7
	EQUIPAGE, #J7
	EQUIPAGE, #C7
	#EQUIPAGE, #C7
	PLATINE, #C3
	PLATINE, #C3
	PLATINE, #C3
	PLATINE, #C3
	CROCHET, #C5
	CROCHET, #C5
	EQUIPAGE, #C6
	EQUIPAGE, #C6
	EQUIPAGE, #C6
	EQUIPAGE, #C6
	##PROBE, #FC#
	#BARRE, #FC#
	#BARRE, #FC#
	#BARRE, #FC#
	#BARRE, #FC#
	#EQUIPAGE, #FC#
	#BARRE, #FC#
	#BARRE, #FC#
	#EQUIPAGE, #S71
	#EQUIPAGE, #S72
	##BIELLETTE,#S6
	##BIELLETTE,#S6
	#PLATINE, #S11
	#BARRE,#S12
	#MANETTE, #S1
	#CROCHET, #S4
	#CROCHET, #S4
	#BARRE,#S5
	#BARRE,#S5
	#PLATINE, #S8
	#PLATINE, #S8
	#CROCHET, #Sx
	BARRE,
	BARRE,
	BARRE, #FC#
	PROBE


])

######NAME OF CONTACT BODY1  DESCRIPTIONS#################################################################
contactBody2=array.array('i',[
	BATIE,#J1
	BATIE,#J1
	BATIE,#C1
	BATIE,#C1
	BATIE,#C1
	BATIE,#C1
	BATIE,#S3
	BATIE,#S3
	AIGUILLE,
	AIGUILLE,
	##BATIE,#J2
	##BATIE,#J2
	PLATINE,#J3
	PLATINE,#J3
	PLATINE,#J4
	PLATINE,#J4
	BARRE,#C4
	BARRE,#C4
	BARRE,#C4
	BARRE,#C4
	BIELLETTE, #J5
	BIELLETTE, #J5
	BIELLETTE, #J6
	BIELLETTE, #J6
	AIGUILLE,
	AIGUILLE,
	##BATIE,#J7
	##BATIE,#J7
	BATIE,#C7
	#BATIE,#C7
	BATIE,#C3
	BATIE,#C3
	BATIE,#C3
	BATIE,#C3
	BARRE,#C5
	BARRE,#C5
	PLATINE,#C6
	PLATINE,#C6
	PLATINE,#C6
	PLATINE,#C6
	##BATIE, #FC#
	#BATIE, #FC#
	#BATIE, #FC#
	#BATIE, #FC#
	#BATIE, #FC#
	#BATIE, #FC#	
	#BATIE, #FC#
	#BATIE, #FC#
	#BATIE, #S71
	#BATIE, #S72
	##BATIE,#S6
	##BATIE,#S6
	#BATIE,#S11
	#BATIE,#S12
	#BATIE,#S1
	#PLATINE,#S4
	#PLATINE,#S4
	#PLATINE,#S5
	#PLATINE,#S5
	#EQUIPAGE, #S8
	#EQUIPAGE, #S8
	#BATIE,#Sx
	PROBE,
	PROBE,
	BATIE,
	BATIE,

])
######CLEARANCE IN THE JOINTS  DESCRIPTIONS#################################################################
ofset=0.001
clc1=0.02

contactOffset=numpy.array([
	0.5-ofset-clc1,#+0.017-0.05,#J1
	0.5-ofset-clc1,#+0.017-0.05,#J1
	ofset,#C1
	ofset,#C1
	ofset,#C1
	ofset,#C1
	0.3-ofset-0.005,#S3
	0.15-ofset-0.005,#S3
	0.5-ofset-clc1,#+0.017-0.05,#J2
	0.5-ofset-clc1,#+0.017-0.05,#J2
	0.5-ofset-clc1,#+0.022-0.05,#J3
	0.5-ofset-clc1,#+0.022-0.05,#J3
	0.5-ofset-clc1,#+0.01-0.05,#J4
	0.5-ofset-clc1,#3+0.01-0.05,#J4
	0.002,#ofset,#C4
	0.002,#ofset,#C4
	0.002,#ofset,#C4
	0.002,#ofset,#C4
	0.5-ofset-clc1,#+0.009-0.05,#J5
	0.5-ofset-clc1,#+0.009-0.05,#J5
	0.5-ofset-clc1,#+0.015-0.05,#0.2 #J6
	0.5-ofset-clc1,#+0.015-0.05,#J6
	0.5-ofset-clc1,#8+0.017-0.05,#J7
	0.5-ofset-clc1,#8+0.017-0.05,#J7
	0.001,#C7
	#0.001,#C7
	ofset,#C3
	ofset,#C3
	ofset,#C3
	ofset,#C3
	0.01,#C5
	0.01,#C5
	ofset,#C6
	ofset,#C6
	ofset,#C6
	ofset,#C6
	##ofset, #FC#
	#ofset, #FC#
	#ofset, #FC#
	#ofset, #FC#
	#ofset, #FC#
	#ofset, #FC#
	#ofset, #FC#
	#ofset, #FC#
	#ofset,#S71
	#ofset,#S72
	##ofset,#S6
	##ofset,#S6
	#ofset,#S11
	#ofset,#S12
	#ofset,#S1
	#ofset,#S4
	#ofset,#S4
	#ofset,#S5
	#ofset,#S5
	#ofset,#S8
	#ofset,#S8
	#ofset,#Sx
	0.01,0.01,
	ofset,
	ofset,
])
######CONTACT OFFSET BODY1  DESCRIPTIONS#################################################################
contactOffsetP1=array.array('I',[
        1,1, #J1
	1,1,1,1, #C1
        1,1,#S3
        1,1, #J2
	1,1, #J3
	1,1, #J4
        1,1,1,1, #C4
	1,1, #J5
	1,1, #J6
	1,1, #J7
        1,#1, #C7
	1,1,1,1, #C3
	1,1, #C5
	1,1,1,1, #C6
	#1,#1, #FC
	#1,1, #FC
	#1,1, #FC
	#1,1, #FC
	#1,1,1,1,1,1,1,1,1,1,1,1,#1,1
	1,1,1,1
 
])
######CONTACT TYPE DESCRIPTIONS#################################################################
contactType3D=array.array('I',[
        1,1, #J1
	1,1,1,1, #C1
        1,1,#S3
        1,1, #J2
	1,1, #J3
	1,1, #J4
        1,1,1,1, #C4
	1,1, #J5
	1,1, #J6
	1,1, #J7
        1,#1, #C7
	1,1,1,1, #C3
	1,1, #C5
	1,1,1,1, #C6
	#1,#1, #FC
	#1,1, #FC
	#1,1, #FC
	#1,1, #FC
	#1,1,1,1,1,1,1,1,1,1,1,1,#1,1
	1,1,1,1
 
])
######COEFFICIENT OF FRICTION DESCRIPTIONS#################################################################
mu=0.1
mu1=0.2
contactmu=array.array('d',[
        mu,mu,#J1
	mu,mu,mu,mu,#C1
        mu1,mu1,#S3
        mu,mu,#J2
        mu,mu,#J3
        mu,mu,#J4
        mu1,mu1,mu1,mu1,#C4
        mu,mu,#J5
        mu,mu,#J6
        mu,mu,#J7
	mu,#mu5,#C7
	mu,mu,mu,mu,#C3
	mu,mu,#C5
	mu,mu,mu,mu,#C6
	#mu,#mu,#FC
	#mu,mu,#FC
	#mu,mu,#FC
	#mu,mu,#FC
	#mu5,mu5,mu0,mu0,mu0,mu0,mu0,mu0,mu0,mu0,mu0,
	mu,mu,mu,mu
])
######COEFFICIENT OF RESTITUTION  DESCRIPTIONS#################################################################
en1=0.0
contacten=array.array('d',[
        en1,en1,#J1
	en1,en1,en1,en1,#C1
        en1,en1,#S3
        en1,en1,#J2
        en1,en1,#J3
        en1,en1,#J4
        en1,en1,en1,en1,#C4
        en1,en1,#J5
        en1,en1,#J6
        en1,en1,#J7
	en1,en1,#C7
        en1,en1,en1,en1,#C3
        en1,en1,#C5
        en1,en1,en1,en1,#C6
        #en1,en1,en1,#en1,#FC
        #en1,en1,en1,en1,#FC
       #en1,en1,en1,en1,en1,en1,en1,en1,en1,en1,en1,
	en1,en1,en1,en1
])

######SIMULATION PARAMETER DESCRIPTIONS#################################################################
#3D parameters
with3D=1
freqUpdate=100
freqOutput=2
stepNumber=20000
dumpGraphic=1
drawMode=mbtb.MBTB_ARTEFACT_REACTION+mbtb.MBTB_ARTEFACT_NORMAL+mbtb.MBTB_ARTEFACT_P1P2
#Simulation parameters
stepSize=1e-3
withProj=0
withReduced=2
et=0.0
solverTol=1e-8
solverIt=10000

#TSdeactivateYPosThreshold=1e-5
#TSdeactivateYVelThreshold=0.0
#TSactivateYPosThreshold=0.0
#TSactivateYVelThreshold=100

#TSProjectionMaxIteration=200
#TSConstraintTol=1e-06
#TSConstraintTolUnilateral=1e-07
#TSLevelOfProjection=1


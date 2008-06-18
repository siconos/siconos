DESCRIPTION DE L'ARBORESCENCE:
-----------------------------
robot = arm2 ou pa10

robotControl.sci = schéma de controle du robot (arm2 ou pa10)
                   même schéma mais différents paramètres
arm2Control.sxd  = schéma "equivalent" dataflow en open-office 
robotPlot.sci    = fcts de plot du robot 2D et 3D
robotLaunch.sce  = Lancement de la simulation (équivalent à un main)
                   on peut lancer sous scilab par 'exec'
                   ou 'scilab -f robotLaunch.sce'


lib/
     modules.sci = Fonctions de base utilisée pour le controle indépendantes
des robots
     simulation.sci = Fonctions de simulation qui appelle l'intégrateur
     plot.sci    = Fonctions de plot des données de simulation

lib/robot (arm2 ou pa10) Modèle dynamique du robot
 Friction.c, Inertia.c,NLEffects.c, Tags.c, RobotModel.h, Makefile
 => robotmodel.a
 Le programme Maple qui a généré ces fonctions C est dans MapleCodeGeneration/
 ScilabInterface.c et modele.sci font l'interface entre les fcts C et les fcts Scilab

lib/pa10/pa10Jac.sci  = Modele Geométrique et jacobienne du PA10
 
-----------------------------------------------------------------------

POUR FAIRE TOURNER L'EXEMPLE PA10:
----------------------------------
+ Lancer Scilab dans RobotSim>

Dans le fichier, pa10Launch.sce, tu peux changer les CI (durée, q0), la
loi de commande (en commentant/décommentant...)

---> exec pa10Launch.sce

-----------------------------------------------------------------------

CE QU'IL ME RESTE A FAIRE:
---------------------------

Par ordre de priorité:
+ Ecrire la doc
+ Fignoler le code (commentaires, en-tête)
+ Fignoler arm2 (re-génération Tags.c, lien Tags.c, faire tourner)
+ Commande cartésienne arm2





#!/usr/bin/env python
##Copyright 2011-2014 Thomas Paviot (tpaviot@gmail.com)
##
##This file is part of pythonOCC.
##
##pythonOCC is free software: you can redistribute it and/or modify
##it under the terms of the GNU Lesser General Public License as published by
##the Free Software Foundation, either version 3 of the License, or
##(at your option) any later version.
##
##pythonOCC is distributed in the hope that it will be useful,
##but WITHOUT ANY WARRANTY; without even the implied warranty of
##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##GNU Lesser General Public License for more details.
##
##You should have received a copy of the GNU Lesser General Public License
##along with pythonOCC.  If not, see <http://www.gnu.org/licenses/>.
#---------------------Converting interstellar into IOJSON--------------------#
from __future__ import print_function
from OCC import STEPControl
import webbrowser
from OCC.Visualization import Tesselator
import OCC
from time import time
import os
import tempfile


import subprocess
from OCC import VERSION
from OCC.gp import gp_Ax1, gp_Pnt, gp_Dir, gp_Trsf, gp_Quaternion,gp_Vec, gp_XYZ
from OCC.TopLoc import TopLoc_Location
from OCC.Display.SimpleGui import get_backend,init_display
from OCC.STEPControl import STEPControl_Reader, STEPControl_Writer, STEPControl_AsIs
from OCC.Interface import Interface_Static_SetCVal
from OCC.BRep import BRep_Builder
from OCC.BRepBuilderAPI import BRepBuilderAPI_Transform
from OCC.TopoDS import TopoDS_Compound
from OCC.IFSelect import IFSelect_RetDone, IFSelect_ItemsByEntity
import OCC.Graphic3d as Graphic3d
from OCC.Quantity import Quantity_NOC_DARKVIOLET, Quantity_NOC_BLUE1, Quantity_NOC_GREEN, Quantity_NOC_RED, Quantity_NOC_ORANGE, Quantity_NOC_SALMON, Quantity_NOC_YELLOW
import sys
import random
import math
import vtk
from vtk.util import numpy_support
vtkmath = vtk.vtkMath()
from OCC.gp import gp_Ax1, gp_Pnt, gp_Dir, gp_Trsf, gp_Quaternion,gp_Vec, gp_XYZ
from OCC.TopLoc import TopLoc_Location

#from Siconos.Mechanics import IO
import siconos.io.mechanics_io as IO
from Quaternion import Quat
from collections import Counter, defaultdict
from itertools import groupby
from operator import itemgetter
from timeit import timeit

# update_progress() : Displays or updates a console progress bar
## Accepts a float between 0 and 1. Any int will be converted to a float.
## A value under 0 represents a 'halt'.
## A value at 1 or bigger represents 100%
def update_progress(progress):
    barLength = 50 # Modify this to change the length of the progress bar
    status = ""
    if isinstance(progress, int):
        progress = float(progress)
    if not isinstance(progress, float):
        progress = 0
        status = "error: progress var must be float\r\n"
    if progress < 0:
        progress = 0
        status = "Halt...\r\n"
    if progress >= 1:
        progress = 1
        status = "Done...\r\n"
    block = int(round(barLength*progress))
    text = "\rPercent: [{0}] {1}% {2}".format( "#"*block + "-"*(barLength-block), progress*100, status)
    sys.stdout.write(text)
    sys.stdout.flush()



HEADER = """
<head>
    <title>pythonOCC @VERSION@ webgl renderer</title>
    <meta name='Author' content='Thomas Paviot - tpaviot@gmail.com'>
    <meta name='Keywords' content='WebGl,pythonOCC'>
    <meta charset="utf-8">
    <style type="text/css">
        body {
            background-color: @background-color@;
            margin: 0px;
            overflow: hidden;
        }
        #info {
            position: absolute;
            top: 96%;
            width: 96%;
            color: #808080;
            padding: 5px;
            font-family: Monospace;
            font-size: 13px;
            text-align: right;
            opacity: 1;
            }
        #pythonocc_rocks {
            padding: 5px;
            position: absolute;
            left: 1%;
            top: 87%;
            height: 60px;
            width: 305px;
            border-radius: 5px;
            border: 2px solid #f7941e;
            opacity: 0.7;
            font-family: Arial;
            background-color: #414042;
            color: #ffffff;
            font-size: 16px;
            opacity: 0.7;
        }

/*--------------------------Progression Bar (begin)---------------------------*/
        #load {
            width: 0%;
            height: 14px;
            background: url( data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAALCAYAAAC+jufvAAAABGdBTUEAALGPC/xhBQAAAAlwSFlzAAAOwAAADsABataJCQAAABp0RVh0U29mdHdhcmUAUGFpbnQuTkVUIHYzLjUuMTAw9HKhAAAAPklEQVQYV2M48Gvvf4ZDv/b9Z9j7Fcha827Df4alr1b9Z1j4YsV/BuML3v8ZTC/7/GcwuwokrG4DCceH/v8Bs2Ef1StO/o0AAAAASUVORK5CYII=);
            -moz-border-radius: 4px;
            border-radius: 4px;
            position: absolute;
            top: 80%;
            left:2%;
            opacity: 1;
        }

        #noload{
            width: 80%;
            background: url(data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAAAEAAAALCAYAAAC+jufvAAAABGdBTUEAALGPC/xhBQAAAAlwSFlzAAAOwAAADsABataJCQAAABp0RVh0U29mdHdhcmUAUGFpbnQuTkVUIHYzLjUuMTAw9HKhAAAANUlEQVQYVy3EIQ4AQQgEwfn/zwghCMwGh8Tj+8yVKN0d2l00M6i70XsPmdmfu6OIQJmJqooPOu8mqi//WKcAAAAASUVORK5CYII=);
            height: 12px;
            -moz-border-radius: 4px;
            border-radius: 4px;
            border: 1px solid #999999;
            position: absolute;
            top: 80%;
            left: 2%;
            opacity: 0.4;
        }

        #loadText {
            font-family: Consolas;
            font-size: 11px;
            color: #000000;
            position: absolute;
            top: 80%;
            left: 42%;
        }

        #progressTime {
            font-size: 11px;
            position: absolute;
            color: #000000;;
            top: 80%;
            left: 82.5%;
        }

        #buttonBar {
            position: absolute;
            top: 83%;
            left: 2%;
        }

        #inset  {
            width: 200px;
            height: 200px;
            background-color: #fff; /* or transparent; will show through only if renderer alpha: true */
            border: 1px solid black; /* or none; */
            margin: 20px;
            padding: 0px;
            position: absolute;
            right: 0%;
            bottom: 15%;
            z-index: 100;
        }

        #label  {
            width: 90px;
            height: 70px;
            background-color: #fff; /* or transparent; will show through only if renderer alpha: true */
            border: 1px solid black; /* or none; */
            margin: 20px;
            padding: 0px;
            position: absolute;
            right: 0%;
            bottom: 3.8%;
            z-index: 100;
        }
/*---------------------------Progression Bar (end)----------------------------*/

/*-----------------------------Coordinates (begin)----------------------------*/
        #coordinates  {
            width: 87%;
            height: 8%;
            background-color: #fff; /* or transparent; will show through only if renderer alpha: true */
            border: 2px solid green; /* or none; */
            margin: 10px;
            padding: 0px;
            position: absolute;
            left: 1%;
            bottom: 2%;
            z-index: 100;
        }

/*------------------------------Coordinates (end)-----------------------------*/

/*-------------------------------Camera (begin)-------------------------------*/
        #film_camera {
            margin: 20px;
            height: 50px;
            padding: 0px;
            position: absolute;
            right: 13.5%;
            top: 1.5%;
            text-align:center;
        }

        /* the following is related to the camera view menu */
        #menu-demo2, #menu-demo2 ul{
            margin: 20px;
            padding: 0px;
            position: absolute;
            left: 4.2%;
            top: -5.8%;
            text-align:center;
        }

        #menu-demo2 li{
            display:inline-block;
            position:relative;
            border-radius:8px 8px 0 0;
        }

        #menu-demo2 ul li{
            display:inherit;
            border-radius:0;
        }

        #menu-demo2 ul li:hover{
            border-radius:0;
        }

        #menu-demo2 ul li:last-child{
            border-radius:0 0 8px 8px;
        }

        #menu-demo2 ul{
            position:absolute;
            max-height:0;
            left: 0;
            right: 0;
            overflow:hidden;
            -moz-transition: .8s all .3s;
            -webkit-transition: .8s all .3s;
            transition: .8s all .3s;
        }

        #menu-demo2 li:hover ul{
            max-height:23em;
        }

        /* background des liens menus */
        #menu-demo2 li:first-child{
            background-color: #65537A;
            background-image:-webkit-linear-gradient(top, #65537A 0%, #2A2333 100%);
            background-image:linear-gradient(to bottom, #65537A 0%, #2A2333 100%);
        }

        #menu-demo2 li:nth-child(2){
            background-color: #729EBF;
            background-image: -webkit-linear-gradient(top, #729EBF 0%, #333A40 100%);
            background-image:linear-gradient(to bottom, #729EBF 0%, #333A40 100%);
        }

        #menu-demo2 li:nth-child(3){
            background-color: #F6AD1A;
            background-image:-webkit-linear-gradient(top, #F6AD1A 0%, #9F391A 100%);
            background-image:linear-gradient(to bottom, #F6AD1A 0%, #9F391A 100%);
        }

        #menu-demo2 li:nth-child(4){
            background-color: #F6AD1A;
            background-image:-webkit-linear-gradient(top, #F6AD1A 0%, #9F391A 100%);
            background-image:linear-gradient(to bottom, #F6AD1A 0%, #9F391A 100%);
        }

        #menu-demo2 li:last-child{
            background-color: #CFFF6A;
            background-image:-webkit-linear-gradient(top, #CFFF6A 0%, #677F35 100%);
            background-image:linear-gradient(to bottom, #CFFF6A 0%, #677F35 100%);
        }

        /* background des liens sous menus */
        #menu-demo2 li:first-child li{
            background:#2A2333;
        }

        #menu-demo2 li:nth-child(2) li{
            background:#333A40;
        }

        #menu-demo2 li:nth-child(3) li{
            background:#9F391A;
        }

        #menu-demo2 li:nth-child(4) li{
            background:#9F391A;
        }

        #menu-demo2 li:last-child li{
            background:#677F35;
        }

        /* background des liens menus et sous menus au survol */
        #menu-demo2 li:first-child:hover, #menu-demo2 li:first-child li:hover{
            background:#65537A;
        }

        #menu-demo2 li:nth-child(2):hover, #menu-demo2 li:nth-child(2) li:hover{
            background:#729EBF;
        }

        #menu-demo2 li:nth-child(3):hover, #menu-demo2 li:nth-child(3) li:hover{
            background:#F6AD1A;
        }

        #menu-demo2 li:nth-child(4):hover, #menu-demo2 li:nth-child(4) li:hover{
            background:#F6AD1A;
        }

        #menu-demo2 li:last-child:hover, #menu-demo2 li:last-child li:hover{
            background:#CFFF6A;
        }

        /* les a href */
        #menu-demo2 a{
            text-decoration:none;
            display:block;
            padding:8px 45px;
            color:#fff;
            font-family:arial;
        }

        #menu-demo2 ul a{
            padding:8px 0;
        }

        #menu-demo2 li:hover li a{
            color:#fff;
            text-transform:inherit;
        }

        #menu-demo2 li:hover a, #menu-demo2 li li:hover a{
            color:#000;
        }
/*--------------------------------Camera (end)--------------------------------*/
/*-----------------------------Add/remove (begin)-----------------------------*/

        /* the following is related to the camera view menu */
        #menu-demo3, #menu-demo3 ul{
            margin: 20px;
            padding: 0px;
            position: absolute;
            left: 14.4%;
            top: -3.1956%;
            text-align:center;
        }

        #menu-demo3 li{
            display:inline-block;
            position:relative;
            border-radius:8px 8px 0 0;
        }

        #menu-demo3 ul li{
            display:inherit;
            border-radius:0;
        }

        #menu-demo3 ul li:hover{
            border-radius:0;
        }

        #menu-demo3 ul li:last-child{
            border-radius:0 0 8px 8px;
        }

        #menu-demo3 ul{
            position:absolute;
            max-height:0;
            left: 0;
            right: 0;
            overflow:hidden;
            -moz-transition: .8s all .3s;
            -webkit-transition: .8s all .3s;
            transition: .8s all .3s;
        }

        #menu-demo3 li:hover ul{
            max-height:53em;
        }

        /* background des liens menus */
        #menu-demo3 li:first-child{
            background-color: #65537A;
            background-image:-webkit-linear-gradient(top, #65537A 0%, #2A2333 100%);
            background-image:linear-gradient(to bottom, #65537A 0%, #2A2333 100%);
        }

        #menu-demo3 li:nth-child(2){
            background-color: #729EBF;
            background-image: -webkit-linear-gradient(top, #729EBF 0%, #333A40 100%);
            background-image:linear-gradient(to bottom, #729EBF 0%, #333A40 100%);
        }

        #menu-demo3 li:nth-child(3){
            background-color: #F6AD1A;
            background-image:-webkit-linear-gradient(top, #F6AD1A 0%, #9F391A 100%);
            background-image:linear-gradient(to bottom, #F6AD1A 0%, #9F391A 100%);
        }

        #menu-demo3 li:nth-child(4){
            background-color: #F6AD1A;
            background-image:-webkit-linear-gradient(top, #F6AD1A 0%, #9F391A 100%);
            background-image:linear-gradient(to bottom, #F6AD1A 0%, #9F391A 100%);
        }

        #menu-demo3 li:last-child{
            background-color: #CFFF6A;
            background-image:-webkit-linear-gradient(top, #CFFF6A 0%, #677F35 100%);
            background-image:linear-gradient(to bottom, #CFFF6A 0%, #677F35 100%);
        }

        /* background des liens sous menus */
        #menu-demo3 li:first-child li{
            background:#2A2333;
        }

        #menu-demo3 li:nth-child(2) li{
            background:#333A40;
        }

        #menu-demo3 li:nth-child(3) li{
            background:#9F391A;
        }

        #menu-demo3 li:nth-child(4) li{
            background:#9F391A;
        }

        #menu-demo3 li:last-child li{
            background:#677F35;
        }

        /* background des liens menus et sous menus au survol */
        #menu-demo3 li:first-child:hover, #menu-demo3 li:first-child li:hover{
            background:#65537A;
        }

        #menu-demo3 li:nth-child(2):hover, #menu-demo3 li:nth-child(2) li:hover{
            background:#729EBF;
        }

        #menu-demo3 li:nth-child(3):hover, #menu-demo2 li:nth-child(3) li:hover{
            background:#F6AD1A;
        }

        #menu-demo3 li:nth-child(4):hover, #menu-demo3 li:nth-child(4) li:hover{
            background:#F6AD1A;
        }

        #menu-demo3 li:last-child:hover, #menu-demo3 li:last-child li:hover{
            background:#CFFF6A;
        }

        /* les a href */
        #menu-demo3 a{
            text-decoration:none;
            display:block;
            padding:8px 45px;
            color:#fff;
            font-family:arial;
        }

        #menu-demo3 ul a{
            padding:8px 0;
        }

        #menu-demo3 li:hover li a{
            color:#fff;
            text-transform:inherit;
        }

        #menu-demo3 li:hover a, #menu-demo3 li li:hover a{
            color:#000;
        }
/*-------------------------------Add/remove (end)-----------------------------*/

        a {
            color: #f7941e;
            text-decoration: none;
        }

        a:hover {
            color: #ffffff;
        }
    </style>
</head>
"""

BODY = """
<body>
    <div id="container"></div>

    <div id="info">
        WebGL engine by <a href="http://github.com/mrdoob/three.js" target="_blank">three.js</a>
        </div>
    </div>

    <div id="inset"></div>
    <div id="label"><font color="red">X : RED</font> </br> <font color="green">Y : GREEN</font> </br> <font color="blue">Z : BLUE</font></br></div> <!-- to remember which axis is X, which is Y, and so on.. -->

    <table  id="coordinates">
        <tr>
            <td id = "coordinatesListPositionXTest">PositionX: </td>
            <td id = "coordinatesListPositionYTest">PositionY: </td>
            <td id = "coordinatesListPositionZTest">PositionZ: </td>
            <td id = "coordinatesListQuaternionXTest">QuaternionX: </td>
            <td id = "coordinatesListQuaternionYTest">QuaternionY: </td>
            <td id = "coordinatesListQuaternionZTest">QuaternionZ: </td>
            <td id = "coordinatesListQuaternionWTest">QuaternionW: </td>
        </tr>
        <tr>
            <td id = "coordinatesListPositionX"></td>
            <td id = "coordinatesListPositionY"></td>
            <td id = "coordinatesListPositionZ"></td>
            <td id="coordinatesListQuaternionX"></td>
            <td id="coordinatesListQuaternionY"></td>
            <td id="coordinatesListQuaternionZ"></td>
            <td id="coordinatesListQuaternionW"></td>
        </tr>
    </table>

    <ul id="menu-demo3" >
        <li><a id="addAndRemove" href="#">Add and Remove</a>                                                                 <!-- change the camera view -->
    	      <ul>
                <li><a id = "ALotOfCubes" href="#" onclick="addAndRemoveCubes()" >Place cubes</li>    <!-- add/remove black and red cubes in/from the scene: it helps to situate oneself in space -->
                <li><a id = "PointsOfApplication" href="#" onclick="addAndRemovePointsOfApplication()" >Remove Points of Application</li>
                <li><a id = "Arrows" href="#" onclick="addAndRemoveArrows()"> Remove Arrows</li>@MENUOBJECTS@
            </ul>
        </li>
    </ul>
    <ul id="menu-demo2">
        <li><a href="#">Camera</a>                                                                 <!-- change the camera view -->
            <ul>
                <li><a href="#" onclick="cameraview('Ox+')">Ox View : in front</a></li>
                <li><a href="#" onclick="cameraview('Ox-')">Ox View : from behind</a></li>
                <li><a href="#" onclick="cameraview('Oy+')">Oy View : from above</a></li>
                <li><a href="#" onclick="cameraview('Oy-')">Oy View : from below</a></li>
                <li><a href="#" onclick="cameraview('Oz+')">Oz View : from the left</a></li>
                <li><a href="#" onclick="cameraview('Oz-')">Oz View : from the right</a></li>
                <li><a href="#" onclick="cameraview('begin')">Initial View</a></li>
            </ul>
        </li>
    </ul>



<!--------------------------Progression Bar (begin)---------------------------->
    <div id="progressBarControl">
        <div id="load"></div>
        <div id="noload" onclick="clickProgress(this, event)"></div>
        <div id="loadText">Pas de lecture</div>
        <span id="progressTime">00:00</span>

        <div id="buttonBar">
            <img  src="@SHARE_PATH@/img/player_play.png" height=20cm id="player" onclick="clickPlay(this)" />
            <img  src="@SHARE_PATH@/img/player_rewind.png" height=20cm id="playerRewind" onclick="clickRewind()"/>
            <img  src="@SHARE_PATH@/img/player_fastforward.png" height=20cm id="playerFastForward" onclick="clickFastForward()" />
        </div>
    </div>
<!---------------------------Progression Bar (end)----------------------------->

    <script type="text/javascript" src="@SHARE_PATH@/threeJS_libraries/Three.js"></script>
    <script type="text/javascript" src="@SHARE_PATH@/threeJS_libraries/OrbitControls.js"></script>
    <script type="text/javascript" src="@SHARE_PATH@/threeJS_libraries/TrackballControls.js"></script>
    <script type="text/javascript" src="@SHARE_PATH@/threeJS_libraries/JQuery.js"></script>
    <script type="text/javascript" src="@SHARE_PATH@/threeJS_libraries/Stats.js"></script>
    <script src="@REN_PATH@/interstellarVectors.json"></script>
    <script src="@REN_PATH@/interstellar.json"></script>@SCRIPTS@
    @VertexShaderDefinition@
    @FragmentShaderDefinition@
    <script type="text/javascript">


//----------------------Variables Declarations------------------------//
var camera, scene, renderer, stats, loader, objects = {}, container;
var directionalLight, directionalLight1, directionalLight2, directionalLight3, directionalLight4, directionalLight5, directionalLight6;
var targetRotation = 0;
var targetRotationOnMouseDown = 0;
var targetRotationY = 0;
var targetRotationYOnMouseDown = 0;
var mouseX = 0;
var mouseXOnMouseDown = 0;
var mouseY = 0;
var mouseYOnMouseDown = 0;
var moveForward = false;
var moveBackward = false;
var moveLeft = false;
var moveRight = false;
var moveUp = false;
var moveDown = false;
var windowHalfX = window.innerWidth / 2;
var windowHalfY = window.innerHeight / 2;
var interstellar = objectClone(dataInterstellar);                // dataInterstellar is in the interstellar.json file, one can't just do interstellar = dataInterstellar at the risk of obtaining a dynamic copy (like a pointer)
var interstellarInit = objectClone(dataInterstellar);            // static copy of dataInterstellar, this one won't change of size, it's just to rewind
var interstellarVectors = objectClone(dataInterstellarVectors);              // ... dataInterstellarFix (related to not moving objects)....
var interstellarVectorsInit = objectClone(dataInterstellarVectors);          // ... dataInterstellarFix ...
var playingNow = false, rewindingNow = 0, fastForwardingNow = 0;            // indicate if: the animation is on pause/playing, rewinding, fastForwardingNow
var timeStepInit = 10;                                                                                                             // animation's normal timestep
var timeStep = timeStepInit;                                                                                                       // timestep may change if we fast forward, rewind etc...
var container2, camera2, scene2, renderer2, axes2, CANVAS_WIDTH = 200, CANVAS_HEIGHT = 200, CAM_DISTANCE = 300;                    // to follow the frame on the left bottom corner
var raycaster = new THREE.Raycaster(), mouse = new THREE.Vector2(), INTERSECTED, FOLLOWINGAVECTOR = false, VECTORTOFOLLOW;             // interactivity
var listOfObjects =[], listOfListOfObjects = [], listOfTypeOfListOfObjects = [], dictionaryListOfObjects={};                                              // important to get the list of objects to interact with
var arrows = {}, pointsOfApplication = {};                                                                        // arrows in order to represent vectors, points of application...
var lotofCubesRed, lotofCubesBlack, THEREARELOTOFCUBES = false, THEREAREARROWS = true, THEREAREPOINTSOFAPPLICATION = true;                                                                    // add/remove black and red cubes in/from the scene: it helps to situate oneself in space
var firstKeyObjects, firstKeyArrows;
var maxSceneLength;
var setOfKeysObjects, setOfKeysArrows, setOfKeysPointsOfApplication,setOfKeysObjectsInterstellar, setOfKeysArrowsInterstellarVectors, firstKeyObjects, firstKeyArrows, firstKeyPointsOfApplication, firstKeyObjectsInterstellar, firstKeyArrowsInterstellarVectors;                    // in order to get a key related to an object of interstellar (one could add key value pairs which key is a string or something else)







init();
animate();

function init() {
    //-------------------------Dom & Renderer-------------------------//
    container = document.getElementById( 'container' );
    renderer = new THREE.WebGLRenderer();
    renderer.setClearColor( 0x9E9E9E, 1 );
    renderer.setPixelRatio( window.devicePixelRatio );
    renderer.setSize( window.innerWidth, window.innerHeight );
    container.appendChild( renderer.domElement );

    //--------------------------Scene---------------------------------//
    scene = new THREE.Scene();
    scene.add( new THREE.AmbientLight(0x101010));

    //------------Camera, Camera 2, Controls, axes, light and cubes-----------//
    // Must be placed after objects: IMPORTANT!

    //-----------------------Dom 2 & Renderer 2-----------------------//
    container2 = document.getElementById( 'inset' );
    renderer2 = new THREE.WebGLRenderer();
    renderer2.setClearColor( 0xf0f0f0, 1 );
    renderer2.setSize( CANVAS_WIDTH, CANVAS_HEIGHT );
    container2.appendChild( renderer2.domElement );

    //--------------------------Scene 2-------------------------------//
    scene2 = new THREE.Scene();

    //----------------------------Axes 2------------------------------//
    axes2 = new THREE.AxisHelper( 100 );
    scene2.add( axes2 );


    @Uniforms@
    @ShaderMaterialDefinition@

    //-----------------------Fix Objects(begin)-----------------------//@STATICOBJECTS@
    //-----------------------Fix Objects(end)-------------------------//


    //-----------------------Objects(begin)---------------------------//@OBJECTS@
    //-----------------------Objects(end)-----------------------------//


    //--------------------------Arrows (begin)------------------------// @ARROWS@
    //---------------------------Arrows (end)-------------------------//


    //------------------Points of Application(begin)------------------//@POINTSOFAPPLICATION@
    //-------------------Points of Application(end)-------------------//

    //-----------------firstKeyArrows//Objects(begin)-----------------//
    // the user may comment an object for a test, let's say the user takes out the object1: the condition in the requestAnimationFrame
    // CANNOT BE  interstellar[ "1" ][ "positionX" ].length>0 and has to be interstellar[ firstKeyObjects ][ "positionX" ].length>0
    // where firstKeyObjects is the key of an object which is present on the scene. The same applies to the arrows and the points of application.

    setOfKeysObjects = Object.keys( objects );
    if (setOfKeysObjects[0]){
        firstKeyObjects = parseInt(setOfKeysObjects[0]);
        var iobject = 1;
        while( !Number.isInteger(firstKeyObjects) && iobject < setOfKeyObjects.length  ){
            firstKeyObjects = parseInt( setOfKeysObjects[ iobject ] );
            iobject++;
        }
    }

    setOfKeysObjectsInterstellar = Object.keys( interstellar );
    firstKeyObjectsInterstellar = parseInt(setOfKeysObjectsInterstellar[0]);
    var iobjectInterstellar = 1;
    while( !Number.isInteger(firstKeyObjectsInterstellar) && iobjectInterstellar < setOfKeysObjectsInterstellar.length  ){
        firstKeyObjectsInterstellar = parseInt( setOfKeysObjectsInterstellar[ iobjectInterstellar ] );
        iobjectInterstellar++;
    }


    setOfKeysArrows = Object.keys( arrows );
    if (setOfKeysArrows[0]){
        firstKeyArrows = parseInt( setOfKeysArrows[0] );
        var iarrow = 0;
        while( !Number.isInteger(firstKeyArrows) && iarrow < setOfKeyArrows.length ){
            firstKeyArrows = parseInt( setOfKeysArrows[ iarrow ] );
            iarrow++;
        }
    }

    setOfKeysArrowsInterstellarVectors = Object.keys( interstellarVectors );
    if (setOfKeysArrowsInterstellarVectors[0]){
        firstKeyArrowsInterstellarVectors = parseInt( setOfKeysArrowsInterstellarVectors[0] );
        var iarrowInterstellarVectors = 0;
        while( !Number.isInteger(firstKeyArrowsInterstellarVectors) && iarrowInterstellarVectors < setOfKeysArrowsInterstellarVectors.length ){
            firstKeyArrowsInterstellarVectors = parseInt( setOfKeysArrowsInterstellarVectors[ iarrowInterstellarVectors ] );
            iarrowInterstellarVectors++;
        }
    }

    setOfKeysPointsOfApplication = Object.keys( pointsOfApplication );
    if (setOfKeysPointsOfApplication[0]){
        firstKeyPointsOfApplication = parseInt( setOfKeysPointsOfApplication[0] );
        var ipointofapplication = 0;
        while( !Number.isInteger(firstKeyPointsOfApplication) && ipointofapplication < setOfKeysPointsOfApplication.length  ){
            firstKeyPointsOfApplication = parseInt( setOfKeysPointsOfApplication[ ipointofapplication ] );
            ipointofapplication++;
        }
    }



    interstellarLength = interstellarInit[ firstKeyObjectsInterstellar ][ "positionX" ].length;       // the total number of timesteps

    //------------------firstKeyArrows//Objects(end)------------------//

    //--------------------------Camera--------------------------------//
    camera = new THREE.PerspectiveCamera( 70, window.innerWidth / window.innerHeight, 0.1, 1000 );
    cameraview('begin');


    //--------------------------Controls------------------------------//
    controls = new THREE.OrbitControls( camera );


    //--------------------------Camera 2------------------------------//
    camera2 = new THREE.PerspectiveCamera( 50, CANVAS_WIDTH / CANVAS_HEIGHT, 1, 1000 );
    camera2.up = camera.up;    // IMPORTANT!


    //----------------------------Axes--------------------------------//
    var SceneLength = getLengthScene();
    axes = new THREE.AxisHelper( 3 * Math.max( SceneLength.x, SceneLength.y, SceneLength.z ) );
    scene.add( axes );
    console.log( "AxisX : RED \\n AxisY : GREEN \\n AxisZ : BLUE " );
    // AxisX : RED
    // AxisY : GREEN
    // AxisZ : BLUE


    //--------------------------Light---------------------------------//
    initLight(scene); // placing 7 lights at strategic points

    //--------------------------Cubes---------------------------------//
    var SceneLength = getLengthScene();
    SceneLength = Math.max( SceneLength.x, SceneLength.y, SceneLength.z );
    var dimensionCubes = 0.1*SceneLength;
    var stepCubes = pow10floor(SceneLength);
    lotofCubesRed = lotLotOfCubes( dimensionCubes, 0.5*stepCubes, 2, 0xff0000 );
    lotofCubesBlack =  lotLotOfCubes( 1.05*dimensionCubes, stepCubes, 2, 0x000000 );   // IMPORTANT!! to set the size of the black cubes slightly bigger than the red ones so that black cubes turn to green whenthe user interacts with them


    //------------------------Stats-----------------------------------//
    stats = new Stats();
    stats.domElement.style.position = 'absolute';
    stats.domElement.style.top = '0px';
    container.appendChild( stats.domElement );

    //------------------------Interactivity---------------------------//
    document.addEventListener( 'mousemove', onDocumentMouseMove, false );                                       // in order to get the coordinates where the mouse is: used for interactivity


    //--------------------Resize window-------------------------------//
    window.addEventListener( 'resize', onWindowResize, false );                                                 // resize the scene if the window is resized
}

function animate() {
    controls.update();                                                                                          // controls: move around objects of the scene, etc...
    camera2.position.copy( camera.position );                                                                   // camera2 follows the frame of the scene
    camera2.position.sub( controls.target ); // added by @libe
    camera2.position.setLength( CAM_DISTANCE );
    camera2.lookAt( scene2.position );
    setTimeout( function() {                                                                                    // forces the animation to follow the timestep

    	requestAnimationFrame( animate );

        for ( var _id in objects ){                                                                                                                                                // for all moving objects
            if( ( interstellar[ firstKeyObjects ][ "positionX" ].length>0 && playingNow ) || ( interstellar[ firstKeyObjects ][ "positionX" ].length>0 && rewindingNow ) || ( interstellar[ firstKeyObjects ][ "positionX" ].length>0 && fastForwardingNow ) ){                // if animation not paused (playingNow) or animation not finished then we do the following:
                if ( ( Number.isInteger( parseInt(_id) ) ) && ( parseInt(_id) >= 0 ) ) {
                    var quaternionTemp = new THREE.Quaternion();
                    quaternionTemp.set( interstellar[ _id ][ "quaternionX" ][ 0 ], interstellar[ _id ][ "quaternionY" ][ 0 ], interstellar[ _id ][ "quaternionZ" ][ 0 ], interstellar[ _id ][ "quaternionW" ][ 0 ] );                                  // set quaternionTemp to be the quaternion at time t
                    var translationInit1 = new THREE.Vector3( interstellar[ _id ][ "initialPosition" ][ "positionX" ], interstellar[ _id ][ "initialPosition" ][ "positionY" ], interstellar[ _id ][ "initialPosition" ][ "positionZ" ] );                                                       // there is an offset for the position
                    translationInit1.copy( translationOffset(translationInit1, quaternionTemp) );

                    objects[ _id ].quaternion.x = quaternionTemp.x;
                    objects[ _id ].quaternion.y = quaternionTemp.y;
                    objects[ _id ].quaternion.z = quaternionTemp.z;
                    objects[ _id ].quaternion.w = quaternionTemp.w;

                    objects[ _id ].position.x = interstellar[ _id ][ "positionX" ][ 0 ] + translationInit1.x;
                    objects[ _id ].position.y = interstellar[ _id ][ "positionY" ][ 0 ] + translationInit1.y;
                    objects[ _id ].position.z = interstellar[ _id ][ "positionZ" ][ 0 ] + translationInit1.z;
                }
            }
        }

        for ( var _id in interstellar ){                                                                                                                                                // for all moving objects
            if( ( interstellar[ firstKeyObjectsInterstellar ][ "positionX" ].length>0 && playingNow ) || ( interstellar[ firstKeyObjectsInterstellar ][ "positionX" ].length>0 && rewindingNow ) || ( interstellar[ firstKeyObjectsInterstellar ][ "positionX" ].length>0 && fastForwardingNow ) ){                // if animation not paused (playingNow) or animation not finished then we do the following:
                if ( ( Number.isInteger( parseInt(_id) ) ) && ( parseInt(_id) >= 0 ) ) {
                    interstellar[ _id ][ "quaternionX" ].shift();
                    interstellar[ _id ][ "quaternionY" ].shift();
                    interstellar[ _id ][ "quaternionZ" ].shift();
                    interstellar[ _id ][ "quaternionW" ].shift();

                    interstellar[ _id ][ "positionX" ].shift();
                    interstellar[ _id ][ "positionY" ].shift();
                    interstellar[ _id ][ "positionZ" ].shift();

                }
            }
        }

        for ( var _jd in pointsOfApplication ){
            if( ( interstellarVectors[ firstKeyPointsOfApplication ][ "pointOfApplicationX" ].length>0 && playingNow ) || ( interstellarVectors[ firstKeyPointsOfApplication ][ "pointOfApplicationY" ].length>0 && rewindingNow ) || ( interstellarVectors[ firstKeyPointsOfApplication ][ "pointOfApplicationZ" ].length>0 && fastForwardingNow ) ){            // if animation not paused (playingNow) or animation not finished then we do the following:
                if ( Number.isInteger( parseInt(_jd) ) ) {
                    pointsOfApplication[ _jd ].position.set(interstellarVectors[ _jd ]["pointOfApplicationX"][0], interstellarVectors[ _jd ]["pointOfApplicationY"][0], interstellarVectors[ _jd ]["pointOfApplicationZ"][0]);
                }
            }
        }

        for ( var _jd in arrows ){                                                                                                                                                          // for all arrows (representing vectors)
            if( ( interstellarVectors[ firstKeyArrows ][ "pointOfApplicationX" ].length>0 && playingNow ) || ( interstellarVectors[ firstKeyArrows ][ "pointOfApplicationY" ].length>0 && rewindingNow ) || ( interstellarVectors[ firstKeyArrows ][ "pointOfApplicationZ" ].length>0 && fastForwardingNow ) ){            // if animation not paused (playingNow) or animation not finished then we do the following:
                if ( Number.isInteger( parseInt(_jd) ) ) {

                    maxSceneLength = Math.max( maxSceneLength.x, maxSceneLength.y, maxSceneLength.z );                                                                                                                      // the following 6 lines produce a normalized force intensity for the scene depending on the maximum intensity of the j-th vector (on all times), the maxSceneLength
                    var forceStrength = Math.sqrt( Math.pow(interstellarVectors[ _jd ][ "forceDirectionX" ][0],2) + Math.pow(interstellarVectors[ _jd ][ "forceDirectionY" ][0],2) + Math.pow(interstellarVectors[ _jd ][ "forceDirectionZ" ][0],2) );
                    forceStrength /= interstellarVectors[ "absoluteMaximumIntensity" ];
                    var pointY = new THREE.Vector3( interstellarVectors[ _jd ][ "forceDirectionX" ][0], interstellarVectors[ _jd ][ "forceDirectionY" ][0], interstellarVectors[ _jd ][ "forceDirectionZ" ][0] );
                    pointY.normalize();
                    pointY.multiplyScalar( forceStrength );

                    var pointX = new THREE.Vector3( interstellarVectors[ _jd ][ "pointOfApplicationX" ][0], interstellarVectors[ _jd ][ "pointOfApplicationY" ][0], interstellarVectors[ _jd ][ "pointOfApplicationZ" ][0] );                                                              // point of application of the j-th vector, that's to say, the point where the vector's arrow begins
                    pointY.add(pointX);
                    vectorRepresentationUpdate( arrows[ _jd ], pointX, pointY );

                }
            }
        }

        for ( var _jd in interstellarVectors ){                                                                                                                                                          // for all arrows (representing vectors)
            if( ( interstellarVectors[ firstKeyArrowsInterstellarVectors ][ "pointOfApplicationX" ].length>0 && playingNow ) || ( interstellarVectors[ firstKeyArrowsInterstellarVectors ][ "pointOfApplicationY" ].length>0 && rewindingNow ) || ( interstellarVectors[ firstKeyArrowsInterstellarVectors ][ "pointOfApplicationZ" ].length>0 && fastForwardingNow ) ){            // if animation not paused (playingNow) or animation not finished then we do the following:
                if ( Number.isInteger( parseInt(_jd) ) ) {
                    interstellarVectors[ _jd ][ "pointOfApplicationX" ].shift();
                    interstellarVectors[ _jd ][ "pointOfApplicationY" ].shift();
                    interstellarVectors[ _jd ][ "pointOfApplicationZ" ].shift();
                    interstellarVectors[ _jd ][ "forceDirectionX" ].shift();
                    interstellarVectors[ _jd ][ "forceDirectionY" ].shift();
                    interstellarVectors[ _jd ][ "forceDirectionZ" ].shift();

                }
            }
        }




    }, timeStep);
    update();
    render();
    stats.update();
}

function render() {
    @IncrementTime@
    interactive( listOfObjects );
    followingAVector();
    renderer.render( scene, camera );
    renderer2.render( scene2, camera2 );
}

function onWindowResize() {
    camera.aspect = window.innerWidth / window.innerHeight;
    camera.updateProjectionMatrix();
    renderer.setSize( window.innerWidth, window.innerHeight );
}


//--------------Progression bar, buttons, etc.. (begin)---------------//

function objectClone( obj ){                                                     // function objectClone returns a static clone of an object
    var result = ( JSON.parse( JSON.stringify( obj ) ) );
    return result;
}

function arrayClone( arr ) {                                                     // function to make a static clone of a list
    var i, copy;
    if( Array.isArray( arr ) ) {
        copy = arr.slice( 0 );
        for( i = 0; i < copy.length; i++ ) {
            copy[ i ] = arrayClone( copy[ i ] );
        }
        return copy;
    }
    else if( typeof arr === 'object' ) {
        throw 'Cannot clone array containing an object!';
    }
    else {
        return arr;
    }
}

function formatTime(time) {                                                      // function formatTime returns the time left until the end of the animation on a min:secs format
    var hours = Math.floor( time / 3600 );
    var mins  = Math.floor( (time % 3600) / 60 );
    var secs  = Math.floor( time % 60 );

    if ( secs < 10 ) {
        secs = "0" + secs;
    }
    if ( hours ) {
        if ( mins < 10 ) {
            mins = "0" + mins;
        }
        return hours + ":" + mins + ":" + secs; // hh:mm:ss
    }
    else {
        return mins + ":" + secs; // mm:ss
    }
}

function getMousePosition( event ) {
    return {
        x: event.pageX,
        y: event.pageY
    };
}

function getPosition( element ){                                                 // position of an element in the dom are regarding it's parent, function getPosition returns the absolute position of this element in the dom
    var top = 0, left = 0;
    do {
        top  += element.offsetTop;
        left += element.offsetLeft;
    } while ( element = element.offsetParent );
    return { x: left, y: top };
}

function update() {                                                              // function update updates the progress bar depending on the size of the interstellar array
    var duration = interstellarLength;                              // total duration
    var time     = interstellar[ firstKeyObjectsInterstellar ][ "positionX" ].length;         // elapsed time
    var fraction = time / duration;

    if ( rewindingNow ) { var percent  = Math.ceil( fraction * 100 ); }
    else {var percent  = 100 - Math.ceil( fraction * 100 );}

    var load = document.getElementById( 'load' );
    load.style.width = 0.8 * ( percent ) + '%';

    var noload = document.getElementById( 'noload' );
    //noload.style.width = 0.8 * ( 100 - percent ) + '%';
    //noload.style.left = 2 + 0.8 * percent + '%';

    var loadText = document.getElementById( 'loadText' );
    loadText.textContent = Math.ceil( percent ) + '%';

    progressTime = document.getElementById( 'progressTime' )
    progressTime.textContent = formatTime( time );
}

function clickProgress( control, event ) {                                       // function clickProgress changes the interstellar array whenever the clickProgress is clicked somewhere so the user is able to choose the instant of the animation
    var parent = getPosition( control );            // absolute position of the progression bar (here it is the noload bar)
    var target = getMousePosition( event );         // the absolute position of the bar that was clicked
    var x = target.x - parent.x;                    // the position of the bar that was clicked starting at the left side corner of the progression bar
    var wrapperWidth = document.getElementById( 'noload' ).offsetWidth;
    var percent = Math.ceil( ( x / wrapperWidth ) * 100 );
    var indexPlayer = Math.ceil( ( x / wrapperWidth ) * interstellarLength );

    for ( var _id in interstellar ){     // for all moving objects, set the interstellar array to the portion of interstellar that corresponds to where the bar was clicked
        if (    ( Number.isInteger(parseInt(_id)) )     &&  ( parseInt(_id) >= 0 )   ){
            interstellar[ _id ][ "positionX" ] = arrayClone( interstellarInit[ _id ][ "positionX" ].slice( indexPlayer ) );
            interstellar[ _id ][ "positionY" ] = arrayClone( interstellarInit[ _id ][ "positionY" ].slice( indexPlayer ) );
            interstellar[ _id ][ "positionZ" ] = arrayClone( interstellarInit[ _id ][ "positionZ" ].slice( indexPlayer ) );
            interstellar[ _id ][ "quaternionX" ] = arrayClone( interstellarInit[ _id ][ "quaternionX" ].slice( indexPlayer ) );
            interstellar[ _id ][ "quaternionY" ] = arrayClone( interstellarInit[ _id ][ "quaternionY" ].slice( indexPlayer ) );
            interstellar[ _id ][ "quaternionZ" ] = arrayClone( interstellarInit[ _id ][ "quaternionZ" ].slice( indexPlayer ) );
            interstellar[ _id ][ "quaternionW" ] = arrayClone( interstellarInit[ _id ][ "quaternionW" ].slice( indexPlayer ) );
        }
    }

    for ( var _jd in interstellarVectors ){     // for all arrows (representing vectors), set the interstellar array to the portion of interstellar that corresponds to where the bar was clicked
        if ( Number.isInteger( parseInt(_jd) ) ){
            interstellarVectors[ _jd ][ "pointOfApplicationX" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationX" ].slice( indexPlayer ) );
            interstellarVectors[ _jd ][ "pointOfApplicationY" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationY" ].slice( indexPlayer ) );
            interstellarVectors[ _jd ][ "pointOfApplicationZ" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationZ" ].slice( indexPlayer ) );
            interstellarVectors[ _jd ][ "forceDirectionX" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionX" ].slice( indexPlayer ) );
            interstellarVectors[ _jd ][ "forceDirectionY" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionY" ].slice( indexPlayer ) );
            interstellarVectors[ _jd ][ "forceDirectionZ" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionZ" ].slice( indexPlayer ) );
        }
    }

}

function clickPlay( control ){
    if( rewindingNow ){
        var duration = interstellarLength;                              // total duration
        var time     = interstellar[ firstKeyObjectsInterstellar ][ "positionX" ].length;         // elapsed time
        for ( var _id in interstellar ){
            if ( parseInt(_id) && ( parseInt(_id) >= 0 ) ){
                interstellar[ _id ][ "positionX" ] = arrayClone( interstellarInit[ _id ][ "positionX" ].slice( time ) );
                interstellar[ _id ][ "positionY" ] = arrayClone( interstellarInit[ _id ][ "positionY" ].slice( time ) );
                interstellar[ _id ][ "positionZ" ] = arrayClone( interstellarInit[ _id ][ "positionZ" ].slice( time ) );
                interstellar[ _id ][ "quaternionX" ] = arrayClone( interstellarInit[ _id ][ "quaternionX" ].slice( time ) );
                interstellar[ _id ][ "quaternionY" ] = arrayClone( interstellarInit[ _id ][ "quaternionY" ].slice( time ) );
                interstellar[ _id ][ "quaternionZ" ] = arrayClone( interstellarInit[ _id ][ "quaternionZ" ].slice( time ) );
                interstellar[ _id ][ "quaternionW" ] = arrayClone( interstellarInit[ _id ][ "quaternionW" ].slice( time ) );
            }
        }
        for ( var _jd in interstellarVectors ){
            if ( Number.isInteger(parseInt(_jd)) ){
                interstellarVectors[ _jd ][ "pointOfApplicationX" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationX" ].slice( time ) );
                interstellarVectors[ _jd ][ "pointOfApplicationY" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationY" ].slice( time ) );
                interstellarVectors[ _jd ][ "pointOfApplicationZ" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationZ" ].slice( time ) );
                interstellarVectors[ _jd ][ "forceDirectionX" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionX" ].slice( time ) );
                interstellarVectors[ _jd ][ "forceDirectionY" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionY" ].slice( time ) );
                interstellarVectors[ _jd ][ "forceDirectionZ" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionZ" ].slice( time ) );
            }
        }


    }
    fastForwardingNow = 0;
    rewindingNow = 0;
    timeStep = timeStepInit;
    var playerButton = document.getElementById( "player" );
    if ( playingNow ) {
        playerButton.src = "@SHARE_PATH@/img/player_play.png";
        playingNow = false;
    }
    else {
        playerButton.src = "@SHARE_PATH@/img/player_pause.png";
        playingNow = true;
    }
}

function clickRewind(){
    fastForwardingNow = 0;
    if ( rewindingNow===0 ){
        rewindingNow = 1;
        var duration = interstellarLength;    // Duree totale
        var time     = interstellar[ firstKeyObjectsInterstellar ][ "positionX" ].length; // Temps ecoule
        var index    = duration - time;
        for ( var _id in interstellar ){
            if ( ( Number.isInteger( parseInt(_id) ) ) && ( parseInt(_id) >= 0 ) ){
                interstellar[ _id ][ "positionX" ] = arrayClone( interstellarInit[ _id ][ "positionX" ].slice( 0, index ) );
                interstellar[ _id ][ "positionX" ].reverse();
                interstellar[ _id ][ "positionY" ] = arrayClone( interstellarInit[ _id ][ "positionY" ].slice( 0, index ) );
                interstellar[ _id ][ "positionY" ].reverse();
                interstellar[ _id ][ "positionZ" ] = arrayClone( interstellarInit[ _id ][ "positionZ" ].slice( 0, index ) );
                interstellar[ _id ][ "positionZ" ].reverse();
                interstellar[ _id ][ "quaternionX" ] = arrayClone( interstellarInit[ _id ][ "quaternionX" ].slice( 0, index ) );
                interstellar[ _id ][ "quaternionX" ].reverse();
                interstellar[ _id ][ "quaternionY" ] = arrayClone( interstellarInit[ _id ][ "quaternionY" ].slice( 0, index ) );
                interstellar[ _id ][ "quaternionY" ].reverse();
                interstellar[ _id ][ "quaternionZ" ] = arrayClone( interstellarInit[ _id ][ "quaternionZ" ].slice( 0, index ) );
                interstellar[ _id ][ "quaternionZ" ].reverse();
                interstellar[ _id ][ "quaternionW" ] = arrayClone( interstellarInit[ _id ][ "quaternionW" ].slice( 0, index ) );
                interstellar[ _id ][ "quaternionW" ].reverse();
            }
        }
        for ( var _jd in interstellarVectors ) {
            if ( Number.isInteger(parseInt(_jd)) ){
                interstellarVectors[ _jd ][ "pointOfApplicationX" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationX" ].slice( 0, index ) );
                interstellarVectors[ _jd ][ "pointOfApplicationX" ].reverse();
                interstellarVectors[ _jd ][ "pointOfApplicationY" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationY" ].slice( 0, index ) );
                interstellarVectors[ _jd ][ "pointOfApplicationY" ].reverse();
                interstellarVectors[ _jd ][ "pointOfApplicationZ" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationZ" ].slice( 0, index ) );
                interstellarVectors[ _jd ][ "pointOfApplicationZ" ].reverse();
                interstellarVectors[ _jd ][ "forceDirectionX" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionX" ].slice( 0, index ) );
                interstellarVectors[ _jd ][ "forceDirectionX" ].reverse();
                interstellarVectors[ _jd ][ "forceDirectionY" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionY" ].slice( 0, index ) );
                interstellarVectors[ _jd ][ "forceDirectionY" ].reverse();
                interstellarVectors[ _jd ][ "forceDirectionZ" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionZ" ].slice( 0, index ) );
                interstellarVectors[ _jd ][ "forceDirectionZ" ].reverse();
            }
        }
    }
    else if ( rewindingNow === 1 ) {
        rewindingNow = 2;
        timeStep /= 10;
    }
    else if ( rewindingNow === 2 ) {
        rewindingNow = 3;
        timeStep /= 10;
    }
    else {}
}

function clickFastForward(){
    rewindingNow = 0;
    if ( fastForwardingNow === 0 ){
        fastForwardingNow = 1;
        timeStep = timeStepInit;
        var duration = interstellarLength;                              // total duration
        var time     = interstellar[ firstKeyObjectsInterstellar ][ "positionX" ].length;         // elapsed time
        var index    = duration - time;
        for ( var _id in interstellar ){
            if ( ( Number.isInteger( parseInt(_id) ) ) && ( parseInt(_id) >= 0 ) ){
                interstellar[ _id ][ "positionX" ] = arrayClone( interstellarInit[ _id ][ "positionX" ].slice( time ) );
                interstellar[ _id ][ "positionY" ] = arrayClone( interstellarInit[ _id ][ "positionY" ].slice( time ) );
                interstellar[ _id ][ "positionZ" ] = arrayClone( interstellarInit[ _id ][ "positionZ" ].slice( time ) );
                interstellar[ _id ][ "quaternionX" ] = arrayClone( interstellarInit[ _id ][ "quaternionX" ].slice( time ) );
                interstellar[ _id ][ "quaternionY" ] = arrayClone( interstellarInit[ _id ][ "quaternionY" ].slice( time ) );
                interstellar[ _id ][ "quaternionZ" ] = arrayClone( interstellarInit[ _id ][ "quaternionZ" ].slice( time ) );
                interstellar[ _id ][ "quaternionW" ] = arrayClone( interstellarInit[ _id ][ "quaternionW" ].slice( time ) );
            }
        }
        for ( var _jd in interstellarVectors ) {
            if ( Number.isInteger(parseInt(_jd)) ){
                interstellarVectors[ _jd ][ "pointOfApplicationX" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationX" ].slice( time ) );
                interstellarVectors[ _jd ][ "pointOfApplicationY" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationY" ].slice( time ) );
                interstellarVectors[ _jd ][ "pointOfApplicationZ" ] = arrayClone( interstellarVectorsInit[ _jd ][ "pointOfApplicationZ" ].slice( time ) );
                interstellarVectors[ _jd ][ "forceDirectionX" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionX" ].slice( time ) );
                interstellarVectors[ _jd ][ "forceDirectionY" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionY" ].slice( time ) );
                interstellarVectors[ _jd ][ "forceDirectionZ" ] = arrayClone( interstellarVectorsInit[ _jd ][ "forceDirectionZ" ].slice( time ) );
            }
        }
    }
    else if ( fastForwardingNow ===1 ) {
        fastForwardingNow = 2;
        timeStep /= 10;
    }
    else if ( rewindingNow === 2 ){
        fastForwardingNow = 3;
        timeStep /= 10;
    }
    else {}
}
//---------------Progression bar, buttons, etc.. (end)----------------//

//---------------Placing of the camera and lights (begin)-------------//
function getCentroid ( mesh ) {                                                  // function getCentroid returns the position of the center of the mesh in the parent's frame
    mesh.geometry.computeBoundingBox();
    boundingBox = mesh.geometry.boundingBox;
    var x0 = boundingBox.min.x;
    var x1 = boundingBox.max.x;
    var y0 = boundingBox.min.y;
    var y1 = boundingBox.max.y;
    var z0 = boundingBox.min.z;
    var z1 = boundingBox.max.z;

    var bWidth = ( x0 > x1 ) ? x0 - x1 : x1 - x0;
    var bHeight = ( y0 > y1 ) ? y0 - y1 : y1 - y0;
    var bDepth = ( z0 > z1 ) ? z0 - z1 : z1 - z0;

    var centroidX = x0 + ( bWidth / 2 ) + mesh.position.x;
    var centroidY = y0 + ( bHeight / 2 )+ mesh.position.y;
    var centroidZ = z0 + ( bDepth / 2 ) + mesh.position.z;

    return { x : centroidX, y : centroidY, z : centroidZ };
}

function getCentroidAbsolute ( mesh ) {                                          // function getCentroidAbsolute returns the absolute position of the center of the mesh
    var centroidX = 0;
    var centroidY = 0;
    var centroidZ = 0;

    var Mesh = mesh;

    do{
        centroidX += getCentroid(Mesh).x;
        centroidY += getCentroid(Mesh).y;
        centroidZ += getCentroid(Mesh).z;
        Mesh = Mesh.parent;
    } while ( Mesh.parent );

    return { x : centroidX, y : centroidY, z : centroidZ };
}

function getLength( mesh ){                                                      // function getLength returns Length, Width, Height of the mesh
    mesh.geometry.computeBoundingBox();
    boundingBox = mesh.geometry.boundingBox;
    var x0 = boundingBox.min.x;
    var x1 = boundingBox.max.x;
    var y0 = boundingBox.min.y;
    var y1 = boundingBox.max.y;
    var z0 = boundingBox.min.z;
    var z1 = boundingBox.max.z;

    var bWidth = ( x0 > x1 ) ? x0 - x1 : x1 - x0;
    var bHeight = ( y0 > y1 ) ? y0 - y1 : y1 - y0;
    var bDepth = ( z0 > z1 ) ? z0 - z1 : z1 - z0;

    return { x : bWidth, y : bHeight, z : bDepth };
}

function getMinMaxAbsolute( mesh ){                                              // function getMinMaxAbsolute returns the min and the max on x,y,z of the mesh
    var gravicenter = getCentroidAbsolute( mesh );
    var bWidth  = getLength( mesh ).x;
    var bHeight = getLength( mesh ).y;
    var bDepth  = getLength( mesh ).z;

    var x0 = gravicenter.x - 0.5 * bWidth;
    var x1 = gravicenter.x + 0.5 * bWidth;
    var y0 = gravicenter.y - 0.5 * bHeight;
    var y1 = gravicenter.y + 0.5 * bHeight;
    var z0 = gravicenter.z - 0.5 * bDepth;
    var z1 = gravicenter.z + 0.5 * bDepth;

    return { xmin : x0, xmax : x1, ymin : y0, ymax : y1, zmin : z0, zmax : z1 };
}

function getMinMaxAbsoluteScene(){                                               // function getMinMaxAbsoluteScene returns the min and the max on x,y,z of the moving objects of the scene
    var x0 = Infinity;
    var x1 = -Infinity;
    var y0 = Infinity;
    var y1 = -Infinity;
    var z0 = Infinity;
    var z1 = -Infinity;

    for ( var _id in objects ) {
        if ( Number.isInteger( parseInt(_id) ) ){
            if ( inObject( objects[ _id ], scene.children ) ){
                var minMaxAbsoluteObject = getMinMaxAbsolute( objects[ _id ] );
                x0 = ( x0 < minMaxAbsoluteObject.xmin ) ? x0 : minMaxAbsoluteObject.xmin;
                x1 = ( x1 > minMaxAbsoluteObject.xmax ) ? x1 : minMaxAbsoluteObject.xmax;
                y0 = ( y0 < minMaxAbsoluteObject.ymin ) ? y0 : minMaxAbsoluteObject.ymin;
                y1 = ( y1 > minMaxAbsoluteObject.ymax ) ? y1 : minMaxAbsoluteObject.ymax;
                z0 = ( z0 < minMaxAbsoluteObject.zmin ) ? z0 : minMaxAbsoluteObject.zmin;
                z1 = ( z1 > minMaxAbsoluteObject.zmax ) ? z1 : minMaxAbsoluteObject.zmax;
            }
        }
    }

    return { xmin : x0, xmax : x1, ymin : y0, ymax : y1, zmin : z0, zmax : z1 };
}

function getLengthScene(){                                                       // function getLengthScene returns Length, Width, Height of the objects of the scene: the 3 caracteristics lengths of the scene
    var minMaxAbsoluteScene = getMinMaxAbsoluteScene();
    var x0 = minMaxAbsoluteScene.xmin;
    var x1 = minMaxAbsoluteScene.xmax;
    var y0 = minMaxAbsoluteScene.ymin;
    var y1 = minMaxAbsoluteScene.ymax;
    var z0 = minMaxAbsoluteScene.zmin;
    var z1 = minMaxAbsoluteScene.zmax;

    var bWidth  = ( x0 > x1 ) ? x0 - x1 : x1 - x0;
    var bHeight = ( y0 > y1 ) ? y0 - y1 : y1 - y0;
    var bDepth  = ( z0 > z1 ) ? z0 - z1 : z1 - z0;

    return { x : bWidth, y : bHeight, z : bDepth };
}


function cameraview( axisName ){                                                 // function cameraview change the position of the cameraview as one clickes on the cameraview menu
    var SceneLength = getLengthScene();
    var minMaxAbsoluteScene = getMinMaxAbsoluteScene();
    var x0 = minMaxAbsoluteScene.xmin;
    var x1 = minMaxAbsoluteScene.xmax;
    var y0 = minMaxAbsoluteScene.ymin;
    var y1 = minMaxAbsoluteScene.ymax;
    var z0 = minMaxAbsoluteScene.zmin;
    var z1 = minMaxAbsoluteScene.zmax;
    var maxSceneLength = Math.max( SceneLength.x, SceneLength.y, SceneLength.z );

    switch ( axisName ) {
        case "Oy+":                                                 // looking from above Oy
            camera.position.x = 0.5 * ( x0 + x1 );                  // position.x in the middle of the scene
            camera.position.y = y1 + 5 * maxSceneLength;            // position.y above 5 times the caracteristic length of the scene
            camera.position.z = 0.5 * ( z0 + z1 );                  // position.z in the middle of the scene
            camera.lookAt( scene.position );
            break;
        case "Oy-":                                                 // looking from underneath Oy
            camera.position.x = 0.5 * ( x0 + x1 );
            camera.position.y = - ( y1 + 5 * maxSceneLength );
            camera.position.z = 0.5 * ( z0 + z1 );
            camera.lookAt( scene.position );
            break;
        case "Oz+":                                                 // ....................... Oz
            camera.position.x = 0.5 * ( x0 + x1 );
            camera.position.y = 0.5 * ( y0 + y1 );
            camera.position.z = z1 + 5 * maxSceneLength;
            camera.lookAt( scene.position );
            break;
        case "Oz-":                                                 // ....................... Oz
            camera.position.x = 0.5 * ( x0 + x1 );
            camera.position.y = 0.5 * ( y0 + y1 );
            camera.position.z = -( z1 + 5 * maxSceneLength );
            camera.lookAt( scene.position );
            break;
        case "Ox+":                                                 // ....................... Ox
            camera.position.x = x1 + 5 * maxSceneLength;
            camera.position.y = 0.5 * ( y0 + y1 );
            camera.position.z = 0.5 * ( z0 + z1 );
            camera.lookAt( scene.position );
            break;
        case "Ox-":                                                 // ....................... Ox
            camera.position.x = -( x1 + 5 * maxSceneLength );
            camera.position.y = 0.5 * ( y0 + y1 );
            camera.position.z = 0.5 * ( z0 + z1 );
            camera.lookAt( scene.position );
            break;
        case "begin":                                               // position of the camera as the animation starts: camera is 4 times the caracteristic length of the scene away from scene
            camera.position.x = x1 + 4 * maxSceneLength;
            camera.position.y = y1 + 4 * maxSceneLength;
            camera.position.z = z1 + 4 * maxSceneLength;
            camera.lookAt( scene.position );
            break;
        default:
            break;
    }
}

function initLight( scene ){                                                     // function initLight add lights to the scene depending on the moving objects of the scene
    var SceneLength = getLengthScene();
    var minMaxAbsoluteScene = getMinMaxAbsoluteScene();
    var x0 = minMaxAbsoluteScene.xmin;
    var x1 = minMaxAbsoluteScene.xmax;
    var y0 = minMaxAbsoluteScene.ymin;
    var y1 = minMaxAbsoluteScene.ymax;
    var z0 = minMaxAbsoluteScene.zmin;
    var z1 = minMaxAbsoluteScene.zmax;
    var maxSceneLength = Math.max( SceneLength.x, SceneLength.y, SceneLength.z );

    directionalLight = new THREE.DirectionalLight( 0xffffff );        // position of the directionalLight as the animation starts: camera is 6 times the caracteristic length of the scene away from scene
    directionalLight.position.x = x1 + 6 * maxSceneLength;
    directionalLight.position.y = y1 + 6 * maxSceneLength;
    directionalLight.position.z = z1 + 6 * maxSceneLength;
    directionalLight.position.normalize();
    scene.add( directionalLight );
// the next source light are placed on the axis Ox, Oy, Oz 3 times the caracteristic length of the scene away from scene

    directionalLight1 = new THREE.DirectionalLight( 0xffffff );
    directionalLight1.position.x = x0 - 3 * maxSceneLength;
    directionalLight1.position.y = 0;
    directionalLight1.position.z = 0;
    directionalLight1.position.normalize();
    scene.add( directionalLight1 );

    directionalLight2 = new THREE.DirectionalLight( 0xffffff );
    directionalLight2.position.x = x1 + 3 * maxSceneLength;
    directionalLight2.position.y = 0;
    directionalLight2.position.z = 0;
    directionalLight2.position.normalize();
    scene.add( directionalLight2 );

    directionalLight3 = new THREE.DirectionalLight( 0xffffff );
    directionalLight3.position.x = 0;
    directionalLight3.position.y = y0 - 3 * maxSceneLength;
    directionalLight3.position.z = 0;
    directionalLight3.position.normalize();
    scene.add( directionalLight3 );

    directionalLight4 = new THREE.DirectionalLight( 0xffffff );
    directionalLight4.position.x = 0;
    directionalLight4.position.y = y1 + 3 * maxSceneLength;
    directionalLight4.position.z = 0;
    directionalLight4.position.normalize();
    scene.add( directionalLight4 );

    directionalLight5 = new THREE.DirectionalLight( 0xffffff );
    directionalLight5.position.x = 0;
    directionalLight5.position.y = 0;
    directionalLight5.position.z = z0 - 3 * maxSceneLength;
    directionalLight5.position.normalize();
    scene.add( directionalLight5 );

    directionalLight6 = new THREE.DirectionalLight( 0xffffff );
    directionalLight6.position.x = 0;
    directionalLight6.position.y = 0;
    directionalLight6.position.z = z1 + 3 * maxSceneLength;
    directionalLight6.position.normalize();
    scene.add( directionalLight6 );
}
//----------------Placing of the camera and lights (end)--------------//

//--------------------Representing vectors (begin)--------------------//

function cylinderMesh( pointX, pointY, radiusBottom, radiusTop, material ) {     // given a pointX of start, a pointY of end, other parameters: radius at the bottom, radius at the top, the material to use, returns a mesh cylinder (but does not add it to the scene)
    var direction = new THREE.Vector3().subVectors( pointY, pointX );
    var orientation = new THREE.Matrix4();
    orientation.lookAt( pointX, pointY, new THREE.Object3D().up );
    orientation.multiply( new THREE.Matrix4().set( 1, 0, 0, 0,
        0, 0, 1, 0,
        0, -1, 0, 0,
        0, 0, 0, 1) );

    var edgeGeometry = new THREE.CylinderGeometry(radiusBottom, radiusTop, direction.length(), 32, 1);
    var edge = new THREE.Mesh(edgeGeometry, material);
    edge.applyMatrix(orientation);
    // position based on midpoints - there may be a better solution than this
    edge.position.x = (pointY.x + pointX.x) / 2;
    edge.position.y = (pointY.y + pointX.y) / 2;
    edge.position.z = (pointY.z + pointX.z) / 2;

    return edge;
}

function vectorRepresentation( pointX, pointY, radiusBottom, radiusTop1, radiusTop2, arrowSmallHeigth , material ) {                // given a pointX of start, a pointY of end, other parameters, returns a mesh arrow (supposed to represent a vector) (but does not add it to the scene)
    var direction = new THREE.Vector3().subVectors( pointY, pointX );
    var orientation = new THREE.Matrix4();
    orientation.lookAt( pointX, pointY, new THREE.Object3D().up );
    var orienInit = new THREE.Vector3().subVectors( pointY, pointX );
    orientation.multiply( new THREE.Matrix4().set( 1, 0, 0, 0,
        0, 0, 1, 0,
        0, -1, 0, 0,
        0, 0, 0, 1) );

    var group = new THREE.Group();                                         // the group will contain the cylinder and the small cone which gathered together form an arrow
    var quatInit =new THREE.Quaternion().copy(group.quaternion);

    var arrowLargeGeometry = new THREE.CylinderGeometry(radiusBottom, radiusBottom, direction.length(), 32, 1);
    var arrowLarge = new THREE.Mesh(arrowLargeGeometry, material);
    arrowLarge.typeOfObjectForInteractivity = "arrow";
    listOfObjects.push(arrowLarge);
    group.add(arrowLarge);

    var arrowSmallGeometry = new THREE.CylinderGeometry(radiusTop1, radiusTop2, arrowSmallHeigth, 32, 1);
    var arrowSmall = new THREE.Mesh(arrowSmallGeometry, material);
    arrowSmall.position.y += direction.length()/2 + arrowSmallHeigth/2;     // placing the small cone at the top of the cylinder, before orientating the arrow, the vector is oriented on the Oy direction
    arrowSmall.typeOfObjectForInteractivity = "arrow"
    listOfObjects.push(arrowSmall);
    group.add(arrowSmall);

    group.applyMatrix(orientation);                                         // apply matrix of orientation to the center of the arrow
    group.position.x = (pointY.x + pointX.x) / 2;                           // place the center of the arrow between pointX and pointY
    group.position.y = (pointY.y + pointX.y) / 2;
    group.position.z = (pointY.z + pointX.z) / 2;

    group.directionToWhere = { x : direction.x, y : direction.y, z : direction.z };  // metadata to the vector so that when mouse on the arrow, it's essential information (point of application, direction) is returned on the screen
    group.fromWhere = { x : pointX.x, y : pointX.y, z : pointX.z };

    var result = { "arrowGroup" : group,
                   "length" : direction.length(),
                   "initialQuaternion" : quatInit,
                   "orientation": orientation,
                   "arrowSmallHeigth" : arrowSmallHeigth };
    return result;        // does not just return the group (representing the arrow) but also other required datas to update the arrow
}

function vectorRepresentationUpdate( arrowMeshDict, pointX, pointY ){           // given a pointX of start, a pointY of end, an existing arrow mesh, updates the orientation and the position of the arrow mesh
    var direction = new THREE.Vector3().subVectors( pointY, pointX );
    var arrowMeshLength = arrowMeshDict[ "length" ];

    arrowMesh = arrowMeshDict[ "arrowGroup" ];
    arrowMesh.quaternion.set( arrowMeshDict[ "initialQuaternion" ].x, arrowMeshDict[ "initialQuaternion" ].y, arrowMeshDict[ "initialQuaternion" ].z, arrowMeshDict[ "initialQuaternion" ].w );
    arrowMesh.position.set( 0, 0, 0 );

    arrowMesh.directionToWhere = { x : direction.x, y : direction.y, z : direction.z };         // updates the metadata
    arrowMesh.fromWhere = { x : pointX.x, y : pointX.y, z : pointX.z };

    var orientation = arrowMeshDict[ "orientation" ];
    orientation.getInverse( orientation );                            // before using an orientation matrix, important to re-initialize on the Oy direction the vector mesh
    arrowMesh.applyMatrix( orientation );
    orientation.lookAt( pointX, pointY, new THREE.Object3D().up );

    orientation.multiply( new THREE.Matrix4().set( 1, 0, 0, 0,
        0, 0, 1, 0,
        0, -1, 0, 0,
        0, 0, 0, 1) );
    arrowMesh.applyMatrix( orientation );

    arrowMesh.position.x = ( pointX.x + pointY.x ) / 2 ;
    arrowMesh.position.y = ( pointX.y + pointY.y ) / 2 ;
    arrowMesh.position.z = ( pointX.z + pointY.z ) / 2 ;

// now the arrow is well oriented and the center of the vector is in middle of pointX and pointY BUT!! the objective is to place the base of the cylinder at the pointX: there are 2 distinct cases   distance between pointX and pointY > arrow's height

    if( direction.length() > arrowMeshLength ){
        pointM = new THREE.Vector3( arrowMesh.position.x, arrowMesh.position.y, arrowMesh.position.z );   // pointM at the center of the arrow
        direction.normalize();                                                                  // normalized vector to go from pointX to pointY, that's to say, at this moment of the code, the normalized vector to go from the bottom of the arrow to the top of the arrow
        direction.multiplyScalar( 0.5 * arrowMeshLength );
        pointM.sub( direction );                                                                // pointM at the bottom of the arrow
        var offset = pointM.distanceTo( pointX );                                               // distance between the bottom of the arrow (pointM) and the pointX
        direction.normalize;
        direction.multiplyScalar( offset );
        arrowMesh.position.x -= 2 * direction.x;
        arrowMesh.position.y -= 2 * direction.y;
        arrowMesh.position.z -= 2 * direction.z;
    }
    else{
        var pointM = [ arrowMesh.position.x, arrowMesh.position.y, arrowMesh.position.z ];
        pointM = new THREE.Vector3( pointM[0], pointM[1], pointM[2] );
        direction.normalize();
        direction.multiplyScalar( 0.5 * arrowMeshLength );
        pointM.sub( direction );
        var offset = pointM.distanceTo( pointX );
        direction.normalize;
        direction.multiplyScalar( offset );
        arrowMesh.position.x +=  2*direction.x;
        arrowMesh.position.y +=  2*direction.y;
        arrowMesh.position.z +=  2*direction.z;
    }
}

function pow10floor(x){                                                          // given a number x, function pow10floor returns the closest power of 10
    return Math.pow( 10, Math.floor( Math.log10(x)+0.5) );
}

function lotLotOfCubes( sizeOfCubes, number, border, coulor ){                   // function lotLotOfCubes creates a list of cubes mesh depending of the network, dimensions to set the cubes, ... (but does not add those cubes to the scene)
    var lotLotOfCubes = [] ;

    var i = -border;
    while ( i <= border ) {
        var j = -border ;
        while ( j <= border) {
            var k = -border;
            while ( k <= border) {
                var cube = new THREE.Mesh( new THREE.BoxGeometry( sizeOfCubes, sizeOfCubes, sizeOfCubes ), new THREE.MeshLambertMaterial({ color : coulor }) );
                cube.castShadow = true;
                cube.position.x = number * i;
                cube.position.y = number * j;
                cube.position.z = number * k;
                cube.typeOfObjectForInteractivity = "cube";
                listOfObjects.push( cube );
                lotLotOfCubes.push( cube );
                k++;
            }
            j++;
        }
        i++;
    }

    return lotLotOfCubes;
}

function loadingToSceneFromArray( arrayOfShapes ){                                       // given a list of meshs, function loadingToSceneFromArray add to the scene all the meshs present in it
    if ( Array.isArray( arrayOfShapes ) ){
        var arrayLength = arrayOfShapes.length;
        for (var iter = 0; iter < arrayLength; iter++){
            scene.add( arrayOfShapes[ iter ] );
        }
    }
    else {
        for (var iter in arrayOfShapes){
            scene.add( arrayOfShapes[ iter ] );
        }
    }
}

function takeOutOfSceneFromArray( arrayOfShapes ){                               // given a list of meshs, function loadingToSceneFromArray remove from the scene all the meshs present in it
    if ( Array.isArray( arrayOfShapes ) ){
        var arrayLength = arrayOfShapes.length;
        for ( var iter = 0; iter < arrayLength; iter++ ) {
            scene.remove( arrayOfShapes[ iter ] );
        }
    }
    else {
        for (var iter in arrayOfShapes){
            scene.remove( arrayOfShapes[ iter ] );
        }
    }

}

function addAndRemovePointsOfApplication(){
    if (THEREAREPOINTSOFAPPLICATION) {
        takeOutOfSceneFromArray(pointsOfApplication);
        var button = document.getElementById("PointsOfApplication");
        button.innerHTML = "Replace Points Of Application";
        THEREAREPOINTSOFAPPLICATION = false;
    }
    else{
        loadingToSceneFromArray(pointsOfApplication);
        var button = document.getElementById("PointsOfApplication");
        button.innerHTML = "Remove Points Of Application";
        THEREAREPOINTSOFAPPLICATION = true;
    }
}

function addAndRemoveArrows(){                                                    // function addAndRemoveArrows add and remove the arrows representing vectors from the scene as the user clicks on the remove/add button
    if (THEREAREARROWS) {                                   // if there are cubes than the user clicked on "remove the arrows"
        for (var iter in arrows){
            scene.remove( arrows[ iter ][ "arrowGroup" ] );
        }                                               // remove from the scene all the arrows from the scene
        var button = document.getElementById( "Arrows" );    // assign variable button to the button responsible to add/remove arrow from the scene
        button.innerHTML = "Replace arrows";                // change the text on the button
        THEREAREARROWS = false;                             // now, there is no more arrow anymore
    }
    else{
        for (var iter in arrows){
            scene.add( arrows[ iter ][ "arrowGroup" ] );
        }
        var button = document.getElementById("Arrows");
        button.innerHTML = "Remove arrows";
        THEREAREARROWS = true;
    }
}

function inObject( element, object ){  // check if element is present in object
    if ( Array.isArray(object) ){
        var count = object.length;
        for ( var iter = 0; iter < count; iter++ ){
            if( object[ iter ] === element ){ return true; }
        }
        return false;
    }
    else{
        for ( var iter in object ){
            if( object[ iter ] === element ){ return true; }
        }
        return false;
    }

}

function addAndRemoveObject( _id ){
    if ( inObject(objects[ _id ], scene.children) ) {
        scene.remove( objects[ _id ] );
        var button = document.getElementById( "Object"+_id.toString() );
        button.innerHTML = "Replace Object "+_id.toString();                    // now, there is no more cubes anymore
    }
    else {
        scene.add( objects[ _id ] );
        var button = document.getElementById( "Object"+_id.toString() );
        button.innerHTML = "Remove Object "+_id.toString();
    }
}

function addAndRemoveCubes(){                                                    // function addAndRemoveCubes add and remove the cubes from the scene as the user clicks on the remove/add button
    if (THEREARELOTOFCUBES) {                                   // if there are cubes than the user clicked on "remove the cubes"
        takeOutOfSceneFromArray(lotofCubesRed);                 // remove from the scene all the red cubes from the scene
        takeOutOfSceneFromArray(lotofCubesBlack);               // ..............................black...................
        var button = document.getElementById("ALotOfCubes");    // assign variable button to the button responsible to add/remove cubes from the scene
        button.innerHTML = "Replace cubes";                // change the text on the button
        THEREARELOTOFCUBES = false;                             // now, there is no more cubes anymore
    }
    else{
        loadingToSceneFromArray(lotofCubesRed);
        loadingToSceneFromArray(lotofCubesBlack);
        var button = document.getElementById("ALotOfCubes");
        button.innerHTML = "Remove cubes";
        THEREARELOTOFCUBES = true;
    }
}

function conjugateQuaternion( quaternionTemp ){                                  // given a quaternion, function conjugateQuaternion returns its conjugate
    var quaternionConjugated = new THREE.Quaternion().copy( quaternionTemp );
    quaternionConjugated.x = (-1)* quaternionConjugated.x;
    quaternionConjugated.y = (-1)* quaternionConjugated.y;
    quaternionConjugated.z = (-1)* quaternionConjugated.z;
    quaternionConjugated.w = quaternionConjugated.w;

    return quaternionConjugated;
}

function translationOffset( offsetTranslation, quaternionTemp ){                 // function translationOffset is a function important to the translation offset
/* given a vector of position [px, py, pz] and a quaternion q = [qw, qx, qy, qz] makes the quaternion multiplication q*[0, px, py, pz]*conj(q)
and returns the vector taken from the x,y,z part of the resulting quaternion multiplication */

    var translationInit1Quaternion = new THREE.Quaternion();
    translationInit1Quaternion.x = offsetTranslation.x;
    translationInit1Quaternion.y = offsetTranslation.y;
    translationInit1Quaternion.z = offsetTranslation.z;
    translationInit1Quaternion.w = 0;


    var quaternionTempConjugated = conjugateQuaternion(quaternionTemp);
    var quaternionResult = new THREE.Quaternion().copy( quaternionTemp );
    quaternionResult.multiplyQuaternions ( quaternionResult, translationInit1Quaternion );
    quaternionResult.multiplyQuaternions( quaternionResult, quaternionTempConjugated );


    var vectorResult = new THREE.Vector3();
    vectorResult.x = quaternionResult.x;
    vectorResult.y = quaternionResult.y;
    vectorResult.z = quaternionResult.z;
    return vectorResult;
}
//---------------------Representing vectors (end)---------------------//



//------------------------Interactivity (begin)-----------------------//

function onDocumentMouseMove( event ) {
	event.preventDefault();
	mouse.x = ( event.clientX / window.innerWidth ) * 2 - 1;
	mouse.y = - ( event.clientY / window.innerHeight ) * 2 + 1;
}


function followingAVector(){
// if mouse overflies a vector on the scene, than it's coordinates will be displayed in real time on the left of the window, until another arrow is overflied

    if(FOLLOWINGAVECTOR){
        var coordinates = document.getElementById( "vectorCoordinatesListDirectionX" );
        coordinates.innerText = VECTORTOFOLLOW.directionToWhere.x;

        coordinates = document.getElementById( "vectorCoordinatesListDirectionY" );
        coordinates.innerText = VECTORTOFOLLOW.directionToWhere.y;

        coordinates = document.getElementById( "vectorCoordinatesListDirectionZ" );
        coordinates.innerText = VECTORTOFOLLOW.directionToWhere.z;

        coordinates = document.getElementById( "vectorCoordinatesListFromWhereX" );
        coordinates.innerText = VECTORTOFOLLOW.fromWhere.x;

        coordinates = document.getElementById( "vectorCoordinatesListFromWhereY" );
        coordinates.innerText = VECTORTOFOLLOW.fromWhere.y;

        coordinates = document.getElementById( "vectorCoordinatesListFromWhereZ" );
        coordinates.innerText = VECTORTOFOLLOW.fromWhere.z;
    }
}

function stopFollowingVector(){
// if the user clicks on the stopFollowingVector button than the current arrow vector is not followed anymore

    FOLLOWINGAVECTOR = false;
    var coordinates = document.getElementById( "vectorCoordinatesListDirectionX" );
    coordinates.innerText = "";

    coordinates = document.getElementById( "vectorCoordinatesListDirectionY" );
    coordinates.innerText = "";

    coordinates = document.getElementById( "vectorCoordinatesListDirectionZ" );
    coordinates.innerText = "";

    coordinates = document.getElementById( "vectorCoordinatesListFromWhereX" );
    coordinates.innerText = "";

    coordinates = document.getElementById( "vectorCoordinatesListFromWhereY" );
    coordinates.innerText = "";

    coordinates = document.getElementById( "vectorCoordinatesListFromWhereZ" );
    coordinates.innerText = "";
}

function interactive( listOfObjects ){                                 // function interactive is responsible for the interaction between the objects and the user

    raycaster.setFromCamera( mouse, camera );                         // a ray beam that will intersect objects of the scene at a given position of the mouse
    var intersects = raycaster.intersectObjects( listOfObjects );     // the list of the intersected objects among listOfTypeOfListOfObjects

    if ( intersects.length > 0 ) {                                                              // if there is at least one intersected object
        if ( INTERSECTED != intersects[ 0 ].object ) {                                          // if previous INTERSECTED is not the first intersected object
            if ( INTERSECTED ) INTERSECTED.material.color.setHex( INTERSECTED.currentHex );     // set the previous INTERSECTED color to its currentHex

            INTERSECTED = intersects[ 0 ].object;                                               // INTERSECTED is now the first intersected object
            INTERSECTED.currentHex = INTERSECTED.material.color.getHex();                       // set INTERSECTED.currentHex to the color of INTERSECTED: in order to return the original color of the object after moving the mouse away from the object
            INTERSECTED.material.color.setHex( 0x1B4F08 );                                      // set the INTERSECTED object color to green

            TYPEOFOBJECT = INTERSECTED.typeOfObjectForInteractivity;                            // the "treatment" is different whether the intersected object is a moving object, an array or a cube....
        }

        switch( TYPEOFOBJECT ) {                                                                // depending on the type of the object, the information to be displayed might differ

            case "parent":
                var coordinatesTextHTML = document.getElementById( "coordinatesListPositionXTest" );
                var coordinates = document.getElementById( "coordinatesListPositionX" );
                coordinatesTextHTML.innerHTML = "<strong> Position X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.position.x;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionYTest" );
                coordinates = document.getElementById( "coordinatesListPositionY" );
                coordinatesTextHTML.innerHTML = "<strong> Position Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.position.y;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionZTest" );
                coordinates = document.getElementById( "coordinatesListPositionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Position Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.position.z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionXTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionX" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.quaternion._x;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionYTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionY" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.quaternion._y;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionZTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.quaternion._z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionWTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionW" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion W: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.quaternion._w;
                break;

            case "object":
                var coordinatesTextHTML = document.getElementById( "coordinatesListPositionXTest" );
                var coordinates = document.getElementById( "coordinatesListPositionX" );
                coordinatesTextHTML.innerHTML = "<strong> Position X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.position.x;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionYTest" );
                coordinates = document.getElementById( "coordinatesListPositionY" );
                coordinatesTextHTML.innerHTML = "<strong> Position Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.position.y;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionZTest" );
                coordinates = document.getElementById( "coordinatesListPositionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Position Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.position.z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionXTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionX" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.quaternion._x;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionYTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionY" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.quaternion._y;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionZTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.quaternion._z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionWTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionW" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion W: </strong>";
                coordinates.innerText = intersects[ 0 ].object.quaternion._w;
                break;

            case "cube":
                var coordinatesTextHTML = document.getElementById( "coordinatesListPositionXTest" );
                var coordinates = document.getElementById( "coordinatesListPositionX" );
                coordinatesTextHTML.innerHTML = "<strong> Position X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.position.x;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionYTest" );
                coordinates = document.getElementById( "coordinatesListPositionY" );
                coordinatesTextHTML.innerHTML = "<strong> Position Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.position.y;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionZTest" );
                coordinates = document.getElementById( "coordinatesListPositionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Position Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.position.z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionXTest" );
                coordinates = document.getElementById("coordinatesListQuaternionX");
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionYTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionY" );
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionZTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionZ" );
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionWTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionW" );
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";
                break;

            case "pointOfApplication":
                var coordinatesTextHTML = document.getElementById( "coordinatesListPositionXTest" );
                var coordinates = document.getElementById( "coordinatesListPositionX" );
                coordinatesTextHTML.innerHTML = "<strong> Position X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.position.x;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionYTest" );
                coordinates = document.getElementById( "coordinatesListPositionY" );
                coordinatesTextHTML.innerHTML = "<strong> Position Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.position.y;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionZTest" );
                coordinates = document.getElementById( "coordinatesListPositionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Position Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.position.z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionXTest" );
                coordinates = document.getElementById("coordinatesListQuaternionX");
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionYTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionY" );
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionZTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionZ" );
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionWTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionW" );
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";
                break;

            case "point":
                var coordinatesTextHTML = document.getElementById( "coordinatesListPositionXTest" );
                var coordinates = document.getElementById( "coordinatesListPositionX" );
                coordinatesTextHTML.innerHTML = "<strong> Position X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.point.x;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionYTest" );
                coordinates = document.getElementById( "coordinatesListPositionY" );
                coordinatesTextHTML.innerHTML = "<strong> Position Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.point.y;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionZTest" );
                coordinates = document.getElementById( "coordinatesListPositionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Position Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.point.z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionXTest" );
                coordinates = document.getElementById("coordinatesListQuaternionX");
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionYTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionY" );
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionZTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionZ" );
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionWTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionW" );
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";
                break;

            case "arrow":
                FOLLOWINGAVECTOR = true;
                VECTORTOFOLLOW = INTERSECTED.parent;
                var coordinatesTextHTML = document.getElementById( "coordinatesListPositionXTest" );
                var coordinates = document.getElementById( "coordinatesListPositionX" );
                coordinatesTextHTML.innerHTML = "<strong> Direction X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.directionToWhere.x;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionYTest" );
                coordinates = document.getElementById( "coordinatesListPositionY" );
                coordinatesTextHTML.innerHTML = "<strong> Direction Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.directionToWhere.y;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionZTest" );
                coordinates = document.getElementById( "coordinatesListPositionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Direction Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.directionToWhere.z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionXTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionX" );
                coordinatesTextHTML.innerHTML = "<strong> Point of Application X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.fromWhere.x;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionYTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionY" );
                coordinatesTextHTML.innerHTML = "<strong> Point of Application Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.fromWhere.y;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionZTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Point of Application Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.parent.fromWhere.z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionWTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionW" );
                coordinatesTextHTML.innerHTML = "";
                coordinates.innerText = "";
                break;

            case "firstChild":
                var coordinatesTextHTML = document.getElementById( "coordinatesListPositionXTest" );
                var coordinates = document.getElementById( "coordinatesListPositionX" );
                coordinatesTextHTML.innerHTML = "<strong> Position X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.children[0].position.x;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionYTest" );
                coordinates = document.getElementById( "coordinatesListPositionY" );
                coordinatesTextHTML.innerHTML = "<strong> Position Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.children[0].position.y;

                coordinatesTextHTML = document.getElementById( "coordinatesListPositionZTest" );
                coordinates = document.getElementById( "coordinatesListPositionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Position Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.children[0].position.z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionXTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionX" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion X: </strong>";
                coordinates.innerText = intersects[ 0 ].object.children[0].quaternion._x;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionYTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionY" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion Y: </strong>";
                coordinates.innerText = intersects[ 0 ].object.children[0].quaternion._y;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionZTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionZ" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion Z: </strong>";
                coordinates.innerText = intersects[ 0 ].object.children[ 0].quaternion._z;

                coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionWTest" );
                coordinates = document.getElementById( "coordinatesListQuaternionW" );
                coordinatesTextHTML.innerHTML = "<strong> Quaternion W: </strong>";
                coordinates.innerText = intersects[ 0 ].object.children[0].quaternion._w;
                break;

            default:
                break;
        }
     }
    else {                                                                                  // if there is no object intersected
      if ( INTERSECTED ) INTERSECTED.material.color.setHex( INTERSECTED.currentHex );       // set back the color of the prievous INTERSECTED to its original color
      var coordinatesTextHTML = document.getElementById( "coordinatesListPositionXTest" );
      var coordinates = document.getElementById( "coordinatesListPositionX" );
      coordinatesTextHTML.innerHTML = "";
      coordinates.innerText = "";

      coordinatesTextHTML = document.getElementById( "coordinatesListPositionYTest" );
      coordinates = document.getElementById( "coordinatesListPositionY" );
      coordinatesTextHTML.innerHTML = "";
      coordinates.innerText = "";

      coordinatesTextHTML = document.getElementById( "coordinatesListPositionZTest" );
      coordinates = document.getElementById( "coordinatesListPositionZ" );
      coordinatesTextHTML.innerHTML = "";
      coordinates.innerText = "";

      coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionXTest" );
      coordinates = document.getElementById( "coordinatesListQuaternionX" );
      coordinatesTextHTML.innerHTML = "";
      coordinates.innerText = "";

      coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionYTest" );
      coordinates = document.getElementById( "coordinatesListQuaternionY" );
      coordinatesTextHTML.innerHTML = "";
      coordinates.innerText = "";

      coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionZTest" );
      coordinates = document.getElementById( "coordinatesListQuaternionZ" );
      coordinatesTextHTML.innerHTML = "";
      coordinates.innerText = "";

      coordinatesTextHTML = document.getElementById( "coordinatesListQuaternionWTest" );
      coordinates = document.getElementById( "coordinatesListQuaternionW" );
      coordinatesTextHTML.innerHTML = "";
      coordinates.innerText = "";

      INTERSECTED = null;
    }
}
//-------------------------Interactivity (end)------------------------//


    </script>
</body>
"""

def translateOneItem(L):
    M = [ L[-1] ] + L[:-1]
    return M

def max_occurrences_1a(seq):
    "dict iteritems"
    c = dict()
    for item in seq:
        c[item] = c.get(item, 0) + 1
    return max(c.iteritems(), key=itemgetter(1))

def different_elements(seq):
    c = dict()
    for item in seq:
        c[item] = c.get(item, 0) + 1
    return len(c)

class HTMLHeader(object):
    def __init__(self, background_color='#000000'):
        self._background_color = background_color

    def get_str(self):
        header_str = HEADER.replace('@background-color@', '%s' % self._background_color)
        header_str = header_str.replace('@VERSION@', OCC.VERSION)
        return header_str


class HTMLBody(object):
    def __init__(self, background_color='#000000', vertex_shader=None,
                 fragment_shader=None, uniforms=None, longueur=1, numberOfVectors=0, interstellar={}, numberOfTimeSteps=1, numberStaticObjects = 0, interstellarVectors = {},
                 ren_path=None):
        self._background_color = background_color
        self._vertex_shader = vertex_shader
        self._fragment_shader = fragment_shader
        self._uniforms = uniforms
        self._longueur = longueur
        self._numberOfVectors = numberOfVectors
        self._interstellar = interstellar
        self._numberOfTimeSteps = numberOfTimeSteps
        self._numberStaticObjects = numberStaticObjects
        self._interstellarVectors = interstellarVectors
        self._path = ren_path
    def get_str(self):
        # get the location where pythonocc is running from
        threejs_build_location = os.sep.join([OCC.__path__[0], 'Display', 'WebGl', 'js'])
        body_str = BODY.replace('@Three.jsPath@', '%s' % threejs_build_location)
        body_str = body_str.replace('@background-color@', '%s' % self._background_color)
        body_str = body_str.replace('@VERSION@', OCC.VERSION)
        body_str = body_str.replace('@longueur@', '%d' % self._longueur )
        body_str = body_str.replace('@numberOfVectors@', '%d' % self._numberOfVectors )
        body_str = body_str.replace('@NUMBEROFTIMESTEPS@', '%d' % self._numberOfTimeSteps )
        body_str = body_str.replace('@SHARE_PATH@', '%s' % share_path )
        body_str = body_str.replace('@REN_PATH@', '%s' % self._path )
        # moving objects are added to the scene but their number depends on the HDF5 file, to load objects in the scene the SCRIPTS, ARROWS, ... python variables take strings of code to add on the final JAVASCRIPT code
        SCRIPTS = ""
        ARROWS = ""
        MENUOBJECTS = ""
        POINTSOFAPPLICATION = ""
        STATICOBJECTS = ""
        OBJECTS = ""
        listOfColors = [ "0x3b3b3b", "0x677E52", "0x6B1A6A", "0x495CFF", "0x673300", "0x333b53" ]
        numberOfColors = len(listOfColors)
        for key in (self._interstellarVectors).keys():
            if ( type(key) is int ):
                POINTSOFAPPLICATION += '\n'
                POINTSOFAPPLICATION += '\n    //-----------------------Point of application #'+str( key )+'-----------------------//'
                POINTSOFAPPLICATION += '\n    var materialPointOfApplication = new THREE.MeshBasicMaterial( { color: 0x046380 } );'
                POINTSOFAPPLICATION += '\n    var lengthMax = getLengthScene();'
                POINTSOFAPPLICATION += '\n    var radiusPointOfApplication = (1/160)*Math.max( lengthMax.x, lengthMax.y, lengthMax.z );                  // radius of 1/80 times the caracteristic length'
                POINTSOFAPPLICATION += '\n    var geometryPointOfApplication = new THREE.SphereGeometry( radiusPointOfApplication, 42, 42 );'
                POINTSOFAPPLICATION += '\n    var pointOfApplication = new THREE.Mesh( geometryPointOfApplication, materialPointOfApplication );'
                POINTSOFAPPLICATION += '\n    pointOfApplication.position.set( interstellarVectors['+str( key )+']["pointOfApplicationX"][0], interstellarVectors['+str( key )+']["pointOfApplicationY"][0], interstellarVectors['+str( key )+']["pointOfApplicationZ"][0]);                // interstellarVectors contains the data needed for the representation of the vectors: point of Application, direction, intensity,..'
                POINTSOFAPPLICATION += '\n    pointOfApplication.typeOfObjectForInteractivity = "pointOfApplication";'
                POINTSOFAPPLICATION += '\n    pointsOfApplication[' + str( key ) + '] = pointOfApplication;'
                POINTSOFAPPLICATION += '\n    listOfObjects.push( pointsOfApplication[' + str( key ) + '] );'
                POINTSOFAPPLICATION += '\n    scene.add( pointsOfApplication[' + str( key ) + '] );'
                ARROWS += '\n'
                ARROWS += '\n    //------------------------Arrow #'+str( key )+'------------------------------//'
                ARROWS += '\n    maxSceneLength = getLengthScene();'
                ARROWS += '\n    maxSceneLength = Math.max( maxSceneLength.x, maxSceneLength.y, maxSceneLength.z );'
                ARROWS += '\n    var maxScene = getMinMaxAbsoluteScene();'
                ARROWS += '\n    maxScene = Math.max( maxScene.x, maxScene.y, maxScene.z );'
                ARROWS += '\n'
                ARROWS += '\n    var forceStrength = interstellarVectors[' + str( key ) + '][ "maximumIntensity" ];'
                ARROWS += '\n    forceStrength /= interstellarVectors[ "absoluteMaximumIntensity" ];'
                ARROWS += '\n    forceStrength *= (0.03)*maxSceneLength;'
                ARROWS += '\n'
                ARROWS += '\n    var pointX = new THREE.Vector3( 0, 100*maxScene, 0 );                  // in the begining vectors are not present but their representation cannot change of size so they are placed far far far... away from the scene so they do not seem present'
                ARROWS += '\n    var pointY = new THREE.Vector3( forceStrength, 100*maxScene, 0 );'
                ARROWS += '\n    var cyl_material = new THREE.MeshBasicMaterial( { color: 0x046380 } );'
                ARROWS += '\n    var cylinderHeigth = new THREE.Vector3().subVectors( pointY, pointX ).length();            // height of the cylinder is the norm of pointY - pointX'
                ARROWS += '\n    var cyl_width = 0.1*forceStrength;                    // default line width'
                ARROWS += '\n    var arrow = vectorRepresentation( pointX, pointY, 0.1*cyl_width, 0.1*cyl_width, cyl_width, 0.1*cylinderHeigth , cyl_material );            // returns an arrow ( in a fact a group, may be good to take a look again in the function vectorRepresentation )'
                ARROWS += '\n    arrows[' + str( key ) + '] = arrow;'
                ARROWS += '\n    scene.add( arrows[' + str( key ) + '][ "arrowGroup" ] );'

        for key in (self._interstellar).keys():
            if ( key < 0 ):
                SCRIPTS += '\n    <script src="'+self._path+'/shape_'+str( abs(key) )+'.js"></script>'
                STATICOBJECTS += '\n'
                STATICOBJECTS += '\n    //-----------------------Object (static) #_'+str( abs(key) )+'-----------------------//'
                STATICOBJECTS += '\n    var loader_'+str( abs(key) )+' = new THREE.JSONLoader();'
                STATICOBJECTS += '\n    var object = loader_'+str( abs(key) )+'.parse((new Shape_'+str( abs(key) )+'()).toJSON().data);'
                STATICOBJECTS += '\n    var materialFix = new THREE.MeshPhongMaterial({ transparent: true, opacity: 0.5 });'
                STATICOBJECTS += '\n    object = new THREE.Mesh( object.geometry, materialFix );'
                STATICOBJECTS += '\n    object.position.set( interstellar[ '+str( key )+' ][ "positionX" ], interstellar[ '+str( key )+' ][ "positionY" ], interstellar[ '+str( key )+' ][ "positionZ" ] );'
                STATICOBJECTS += '\n    object.quaternion.set( interstellar[ '+str( key )+' ][ "quaternionX" ], interstellar[ '+str( key )+' ][ "quaternionY" ], interstellar[ '+str( key )+' ][ "quaternionZ" ], interstellar[ '+str( key )+' ][ "quaternionW" ] );'
                STATICOBJECTS += '\n    objects[ '+str(key)+' ] = object;'
                STATICOBJECTS += '\n    scene.add( objects[ '+str(key)+' ] );'

        for key in (self._interstellar).keys():
            if ( key >= 0 ):
                SCRIPTS += '\n    <script src="'+self._path+'/shape'+str(key)+'.js"></script>'
                MENUOBJECTS += '\n              <li><a id = "Object'+str(key)+'" href="#" onclick="addAndRemoveObject('+str(key)+')"> Remove object '+str(key)+'</li>'
                OBJECTS += '\n'
                OBJECTS += '\n    //------------------------Object #'+str(key)+'------------------------------//'
                OBJECTS += '\n    var loader'+str(key)+' = new THREE.JSONLoader();'
                OBJECTS += '\n    var object = loader'+str(key)+'.parse((new Shape'+str(key)+'()).toJSON().data);'
                OBJECTS += '\n    object = new THREE.Mesh( object.geometry, new THREE.MeshPhongMaterial( { transparent: true, opacity: 0.9 } ) );'
                color = listOfColors[0]                                             # we change the colors, a list is permuted everytime a new object is added
                listOfColors = translateOneItem(listOfColors)
                OBJECTS += '\n    object.material.color.setHex('+ color +');'
                OBJECTS += '\n    object.typeOfObjectForInteractivity = "object";'
                OBJECTS += '\n    objects[ '+str(key)+' ] = object;'
                OBJECTS += '\n    objects[ '+str(key)+' ].overdraw = true;'
                #--------  the following code in comment is not needed but the procedure is good and might be good another time, the procedure is to create a cube of dimensions 0, to initially place it on the center of the object, to add the object to the cube and by so, we move the cube of dimensions 0 to move the object  ----------#
                """
                OBJECTS += '\n    /*var gravityCenterX = getCentroid(objects['+str(key)+']).x;'
                OBJECTS += '\n    var gravityCenterY = getCentroid(objects['+str(key)+']).y;'
                OBJECTS += '\n    var gravityCenterZ = getCentroid(objects['+str(key)+']).z;'
                OBJECTS += '\n    var cube = new THREE.Mesh(new THREE.BoxGeometry(0,0,0), new THREE.MeshLambertMaterial());'
                OBJECTS += '\n    cube.position.x = gravityCenterX;'
                OBJECTS += '\n    cube.position.y = gravityCenterY;'
                OBJECTS += '\n    cube.position.z = gravityCenterZ;'
                OBJECTS += '\n    objects[ '+str(key)+' ].position.x -= gravityCenterX;'
                OBJECTS += '\n    objects[ '+str(key)+' ].position.y -= gravityCenterY;'
                OBJECTS += '\n    objects[ '+str(key)+' ].position.z -= gravityCenterZ;'
                OBJECTS += '\n    objects['+str(key)+'].castShadow = true;'
                OBJECTS += '\n    cubes.push( cube );'
                OBJECTS += '\n    cubes[ '+str(key)+' ].add(objects[ '+str(i)+' ]);//*/'
                """
                OBJECTS += '\n'
                OBJECTS += '\n    var quaternionTemp = new THREE.Quaternion();'
                OBJECTS += '\n    quaternionTemp.set( interstellar[ '+str(key)+' ][ "quaternionX" ][ 0 ], interstellar[ '+str(key)+' ][ "quaternionY" ][ 0 ], interstellar[ '+str(key)+' ][ "quaternionZ" ][ 0 ], interstellar[ '+str(key)+' ][ "quaternionW" ][ 0 ] );'
                OBJECTS += '\n    var translationInit1 = new THREE.Vector3( interstellar[ '+str(key)+' ][ "initialPosition" ][ "positionX" ], interstellar[ '+str(key)+' ][ "initialPosition" ][ "positionY" ], interstellar[ '+str(key)+' ][ "initialPosition" ][ "positionZ" ] );'
                OBJECTS += '\n    translationInit1.copy( translationOffset(translationInit1, quaternionTemp) );'
                OBJECTS += '\n'
                OBJECTS += '\n    objects[ '+str(key)+' ].quaternion.x = quaternionTemp.x;'
                OBJECTS += '\n    objects[ '+str(key)+' ].quaternion.y = quaternionTemp.y;'
                OBJECTS += '\n    objects[ '+str(key)+' ].quaternion.z = quaternionTemp.z;'
                OBJECTS += '\n    objects[ '+str(key)+' ].quaternion.w = quaternionTemp.w;'
                OBJECTS += '\n    interstellar[ '+str(key)+' ][ "quaternionX" ].shift();'
                OBJECTS += '\n    interstellar[ '+str(key)+' ][ "quaternionY" ].shift();'
                OBJECTS += '\n    interstellar[ '+str(key)+' ][ "quaternionZ" ].shift();'
                OBJECTS += '\n    interstellar[ '+str(key)+' ][ "quaternionW" ].shift();'
                OBJECTS += '\n'
                OBJECTS += '\n'
                OBJECTS += '\n    objects[ '+str(key)+' ].position.x = interstellar[ '+str(key)+' ][ "positionX" ][ 0 ] + translationInit1.x;'
                OBJECTS += '\n    objects[ '+str(key)+' ].position.y = interstellar[ '+str(key)+' ][ "positionY" ][ 0 ] + translationInit1.y;'
                OBJECTS += '\n    objects[ '+str(key)+' ].position.z = interstellar[ '+str(key)+' ][ "positionZ" ][ 0 ] + translationInit1.z;'
                OBJECTS += '\n    interstellar[ '+str(key)+' ][ "positionX" ].shift();'
                OBJECTS += '\n    interstellar[ '+str(key)+' ][ "positionY" ].shift();'
                OBJECTS += '\n    interstellar[ '+str(key)+' ][ "positionZ" ].shift();'
                OBJECTS += '\n    listOfObjects.push( objects['+str(key)+'] );'
                OBJECTS += '\n    scene.add( objects['+str(key)+'] );'
        body_str = body_str.replace( '@SCRIPTS@', SCRIPTS )
        body_str = body_str.replace( '@ARROWS@', ARROWS )
        body_str = body_str.replace( '@STATICOBJECTS@', STATICOBJECTS )
        body_str = body_str.replace( '@MENUOBJECTS@', MENUOBJECTS )
        body_str = body_str.replace( '@POINTSOFAPPLICATION@', POINTSOFAPPLICATION )
        body_str = body_str.replace( '@OBJECTS@', OBJECTS )


        if (self._fragment_shader is not None) and (self._fragment_shader is not None):
            vertex_shader_string_definition = '<script type="x-shader/x-vertex" id="vertexShader">%s</script>' % self._vertex_shader
            fragment_shader_string_definition = '<script type="x-shader/x-fragment" id="fragmentShader">%s</script>' % self._fragment_shader
            shader_material_definition = """
            var vertexShader = document.getElementById( 'vertexShader' ).textContent;
            var fragmentShader = document.getElementById( 'fragmentShader' ).textContent;
            var shader_material = new THREE.ShaderMaterial( { uniforms: uniforms,
                                                              vertexShader: vertexShader,
                                                              fragmentShader: fragmentShader } );
            """
            if self._uniforms is None:
                body_str = body_str.replace('@Uniforms@', 'uniforms ={};\n')
                body_str = body_str.replace('@IncrementTime@', '')
                #-----------------MODIFICATIONS--------
            else:
                body_str = body_str.replace('@Uniforms@', self._uniforms)
                if 'time' in self._uniforms:
                    backbody_str = body_str.replace('@IncrementTime@', 'uniforms.time.value += 0.05;')
                else:
                    body_str = body_str.replace('@IncrementTime@', '')
            body_str = body_str.replace('@VertexShaderDefinition@', vertex_shader_string_definition)
            body_str = body_str.replace('@FragmentShaderDefinition@', fragment_shader_string_definition)
            body_str = body_str.replace('@ShaderMaterialDefinition@', shader_material_definition)
            body_str = body_str.replace('@ShapeMaterial@', 'shader_material')
        else:
            body_str = body_str.replace('@Uniforms@', '')
            body_str = body_str.replace('@VertexShaderDefinition@', '')
            body_str = body_str.replace('@FragmentShaderDefinition@', '')
            body_str = body_str.replace('@ShaderMaterialDefinition@', '')
            body_str = body_str.replace('@ShapeMaterial@', 'phong_material')
            body_str = body_str.replace('@IncrementTime@', '')
        return body_str


class ThreejsRenderer(object):
    def __init__(self, background_color="#123345", vertex_shader=None, fragment_shader=None, uniforms=None, path=None, longueur=1, numberOfVectors=0,interstellar = {}, numberOfTimeSteps = 1, numberStaticObjects = 0, interstellarVectors={} ):
        if not path:
            self._path =  tempfile.mkdtemp(prefix='.shape')
        else:
            self._path = path

        self._js_filename = os.path.join(self._path, "shape.js")
        self._html_filename = os.path.join(output_path,"index.html"    )
        self._background_color = background_color
        self._vertex_shader = vertex_shader
        self._fragment_shader = fragment_shader
        self._uniforms = uniforms
        self._longueur = longueur
        self._numberOfVectors = numberOfVectors
        self._interstellar = interstellar
        self._numberOfTimeSteps = numberOfTimeSteps
        self._numberStaticObjects = numberStaticObjects
        self._interstellarVectors = interstellarVectors

    def set_vertex_shader(self, vertex_shader):
        ''' adds a vertex shader definition '''
        self._vertex_shader = vertex_shader

    def set_fragment_shader(self, fragment_shader):
        ''' adds a fragment shader '''
        self._fragment_shader = fragment_shader

    def create_files(self, shape, filename=None):
        ''' generate .js and .html files '''
        self._shape = shape
        print("Tesselate shape ...")
        t0 = time()
        tess = Tesselator(self._shape)
        t1 = time()
        print("done in %f s." % (t1-t0))
        print("Exporting tesselation to JSON ...")
        t2 = time()
        if filename == None:
        	tess.ExportShapeToThreejs(self._js_filename)
        else:
            tess.ExportShapeToThreejs(filename)
        t3 = time()
        print("done in %f s." % (t3-t2))
        print("Generating HTML stream ...")
        self.GenerateHTMLFile()
        print("done.")
        return self._js_filename, self._html_filename

    def DisplayShape(self, shape):
        self.create_files(shape)
        print("Opening html output in the default webbrowser ...")
        # previous version us a os.system call to the "open" command
        # but this is a platform (osx) specific solution
        _path = "file:///{0}".format(os.path.join(os.getcwd(), self._html_filename))
        webbrowser.open_new_tab(_path)

    def CreateDictionaryOfShapes(self, dictionaryOfShapes):
        for key in dictionaryOfShapes.keys():
            if ( key >= 0 ):
                filename=os.path.join(self._path, "shape"+str( key )+".js")
                print(key)
            if ( key < 0 ):
                filename=os.path.join(self._path, "shape_"+str( abs(key) )+".js")
            #self._shapeLists.append(filename)
            shape = dictionaryOfShapes[ key ]
            self.create_files(shape,filename)
            print("Creating the html output in the default directory ...")
            # previous version us a os.system call to the "open" command
            # but this is a platform (osx) specific solution
            _path = "file:///{0}".format(os.path.join(os.getcwd(), self._html_filename))
            #os.rename(_path, "/home/inria/siconos/Examples/Mechanics/Mechanisms/SliderCrank/")
            #webbrowser.open_new_tab(_path)

    def GenerateHTMLFile(self):
        """ Generate the HTML file to be rendered wy the web browser
        """
        fp = open(self._html_filename, "w")
        fp.write("<!DOCTYPE HTML>")
        fp.write('<html lang="en">')
        # header
        fp.write(HTMLHeader(self._background_color).get_str())
        # body
        fp.write(HTMLBody(self._background_color,
                          self._vertex_shader,
                          self._fragment_shader,
                          self._uniforms, interstellar = self._interstellar,longueur = self._longueur, numberOfVectors = self._numberOfVectors, numberStaticObjects = self._numberStaticObjects, interstellarVectors = self._interstellarVectors, ren_path = self._path ).get_str())
        fp.write("</html>\n")
        fp.close()
def usage():
    print(' usage :  renderer --help --share_path=xxx --output_path=xxx file.hdf5')

#-----------------------------MODIFICATIONS------------------------------------#
if __name__ == "__main__":

    output_path=os.getcwd()
    bin_path = os.path.dirname(os.path.realpath(__file__))
    share_path = os.path.join(bin_path, "../share/siconos/")
    import sys, getopt
    try:
        if len(sys.argv) < 2:
            raise RuntimeError("an hdf5 file is needed")
        opts, args = getopt.getopt(sys.argv[1:],"",["help","share_path=","output_path="])
        print(opts,args)
    except getopt.GetoptError:
        usage()
        sys.exit(2)
    except RuntimeError as e :
        print("-->", e)
        usage()
        sys.exit(2)
    for opt, arg in opts:
        print(opt,arg)
        if opt == '-h':
            usage()
            sys.exit()
        elif opt in ("--share_path"):
            share_path = arg
        elif opt in ("--output_path"):
            output_path = os.path.join(arg)

    print(sys.argv[-1])
    tempfile.tempdir=output_path
    from OCC.BRepPrimAPI import BRepPrimAPI_MakeBox
    with IO.Hdf5(sys.argv[-1], 'r') as io:    #sys.argv[1]
        nbobjs = len(filter(lambda x: io.instances()[x].attrs['id'] >= 0, io.instances()))                # number of moving objects
        numberStaticObjects = len(filter(lambda x: io.instances()[x].attrs['id'] < 0, io.instances()))    # number of static objects
        longueur = nbobjs
        print("\nNumber of moving objects: %d" % nbobjs)
        print("########### Load shapes")
        obj_by_id = dict()
        #----------------------Loading .step and .stp Files ( begin)------------------------#
        dictionaryOfShapes = {}                                   # dictionaryOfShapes will contain the mesh of the objects
        translation ={}                                     # translation and orientation are dictionaries that will contain the offset of translation and quaternion
        orientation ={}

        print("%s" %dictionaryOfShapes)
        for instance in io.instances():                    # we run through the id's and not the order of the objects so we might see object3 first then object1 and object2 so we need to initialize first the list
            id = io.instances()[instance].attrs['id']
            obj_by_id[id] = instance
            obj = instance

            print("\nobj: %s"  % obj)
            if (id >= 0):
                shape_name  = obj[4:]                    # to get the correct name, it is IMPORTANT to take out a small part of the name IMPORTANT !!!!!
            if (id < 0):
                shape_name  = obj[9:]                    # get the name WE TAKE OUT "ARTEFACT" IMPORTANT!!!!

            print("\nid : %s" % id)
            translation[id] = io.instances()[instance][shape_name+'-0'].attrs['translation']
            orientation[id] = io.instances()[instance][shape_name+'-0'].attrs['orientation']
            print("\ntranslation",repr(translation))
            print("\norientation",repr(orientation))
            print("\nshape_name : %s" % shape_name)
            print("\n%s" % io.instances()[obj])                              # begin of the file: .stp or .step
            with IO.tmpfile(contents=io.shapes()[shape_name][:][0]) as tmpfile:
                step_reader = STEPControl_Reader()
                status = step_reader.ReadFile(tmpfile[1])                          # step FILE
                if status == IFSelect_RetDone:  # check status
                    failsonly = False
                    step_reader.PrintCheckLoad(failsonly, IFSelect_ItemsByEntity)
                    step_reader.PrintCheckTransfer(failsonly, IFSelect_ItemsByEntity)
                step_reader.TransferRoot()
                fichier = step_reader.Shape()
                dictionaryOfShapes[id] = fichier

        print("\n%s" % dictionaryOfShapes)
        print(translation[1])
        #----------------------Loading .step and .stp Files ( end(almost))------------------------#









        #-----------------------------Loading positions----------------------------#
        print("\n ########### Loading positions")
        dpos_data = io.dynamic_data()[:]                  # moving objects
        print("\n%s" % dpos_data)

        print("\n ########### Loading vectors")
        cf_data = io.contact_forces_data()[:]              # vectors
        print("\n%s" % cf_data[0])











        #----------------------Vectors positions orientations----------------------#
        lengthOfTimes = len( cf_data )
        listOfTimes   = []
        for i in range( lengthOfTimes ):
            listOfTimes.append( cf_data[ i ][ 0 ] )
        numberOfVectors = max_occurrences_1a( listOfTimes )[ 1 ]
        numberOfTimeStepsVectors = different_elements( listOfTimes )


        #----------------initialisation of the list containing the vectors
        interstellarVectors = {}
        for vectorNumber in range( numberOfVectors ):
            interstellarVectors[ vectorNumber ] = {}
            interstellarVectors[ vectorNumber ][ "timeStep" ] = []
            interstellarVectors[ vectorNumber ][ "pointOfApplicationX" ] = []
            interstellarVectors[ vectorNumber ][ "pointOfApplicationY" ] = []
            interstellarVectors[ vectorNumber ][ "pointOfApplicationZ" ] = []
            interstellarVectors[ vectorNumber ][ "forceDirectionX" ] = []
            interstellarVectors[ vectorNumber ][ "forceDirectionY" ] = []
            interstellarVectors[ vectorNumber ][ "forceDirectionZ" ] = []

        #---------------We get the positions and the directions of the vectors in a huge list called listOfVectors
        # the important thing here is that the vectors are in file similar to an excel and in this file if the vector is full of 0 then it is not showed that's why the following is used
        i = 0
        rowNumber = 0
        lastRowNumber = 0
        j = 0
        while( (i+j) < numberOfTimeStepsVectors*numberOfVectors ):               # while the whole file is not overflied
            timeStepNow = cf_data[rowNumber][0]                                  # timesteps represents the first arrow at a given time
            j = 0
            while ( timeStepNow == cf_data[rowNumber+j][0] ):                    # while we're still in the same timestep we had the data to listOfVectors
                interstellarVectors[ j ][ "timeStep" ].append( float(cf_data[rowNumber+j][ 0 ]) )
                interstellarVectors[ j ][ "pointOfApplicationX" ].append( float(cf_data[rowNumber+j][ 2 ]) )
                interstellarVectors[ j ][ "pointOfApplicationY" ].append( float(cf_data[rowNumber+j][ 3 ]) )
                interstellarVectors[ j ][ "pointOfApplicationZ" ].append( float(cf_data[rowNumber+j][ 4 ]) )
                interstellarVectors[ j ][ "forceDirectionX" ].append( float(cf_data[rowNumber+j][ 11 ]) )
                interstellarVectors[ j ][ "forceDirectionY" ].append( float(cf_data[rowNumber+j][ 12 ]) )
                interstellarVectors[ j ][ "forceDirectionZ" ].append( float(cf_data[rowNumber+j][ 13 ]) )
                j += 1
            rowNumber += j
            numberOfVectorsLeftToCheck = numberOfVectors - j
            while( numberOfVectorsLeftToCheck > 0 ):                             # during timesteps if we did not overflied all the arrows then we complete by the corresponding numebr of 0 in listOfVectors
                #interstellarVectors[ j ][ "timeStep" ].append( float(cf_data[rowNumber+j][ 0 ]) )
                interstellarVectors[ j ][ "pointOfApplicationX" ].append( 0 )
                interstellarVectors[ j ][ "pointOfApplicationY" ].append( 0 )
                interstellarVectors[ j ][ "pointOfApplicationZ" ].append( 0 )
                interstellarVectors[ j ][ "forceDirectionX" ].append( 0 )
                interstellarVectors[ j ][ "forceDirectionY" ].append( 0 )
                interstellarVectors[ j ][ "forceDirectionZ" ].append( 0 )
                numberOfVectorsLeftToCheck -= 1
                j += 1
            i += numberOfVectors
            lastRowNumber = rowNumber

        # for a given vector arrow it is also important to take the maximum norm of this vector in the entire cycle of the animation
        listOfMaxima = []
        for key in interstellarVectors.keys():
            liste = []
            for i in range( len( interstellarVectors[ key ][ "forceDirectionX" ] ) ):
                u = 0
                u += ( interstellarVectors[ key ][ "forceDirectionX" ][ i ] ) ** 2
                u += ( interstellarVectors[ key ][ "forceDirectionY" ][ i ] ) ** 2
                u += ( interstellarVectors[ key ][ "forceDirectionZ" ][ i ] ) ** 2
                u = math.sqrt(u)
                liste.append(u)
            maxListe = max( liste )
            interstellarVectors[ key ][ "maximumIntensity" ] = maxListe
            listOfMaxima.append( maxListe )
        interstellarVectors["absoluteMaximumIntensity"] = max(listOfMaxima)        # we shall also consider the maximum of the maxima vectors






        #--------------Number of Time Steps, Positions and Rotations---------------#
        numberOfTimeSteps= dpos_data.shape[0]/nbobjs
        interstellar = {}
        interstellarFix = {}
        positions_ini = dpos_data[nbobjs*0:nbobjs*0+nbobjs, 2:]

        for instance in io.instances():
            _id = io.instances()[instance].attrs['id']
            interstellar[_id] = {}
            if (_id >= 0):
                interstellar[_id][ "positionX" ] = []      # position x
                interstellar[_id][ "positionY" ] = []      # position y
                interstellar[_id][ "positionZ" ] = []      # position z
                interstellar[_id][ "quaternionX" ] = []      # quaternion coordoninate x
                interstellar[_id][ "quaternionY" ] = []      # quaternion coordoninate y
                interstellar[_id][ "quaternionZ" ] = []      # quaternion coordoninate z
                interstellar[_id][ "quaternionW" ] = []      # quaternion coordoninate w
                interstellar[_id][ "initialPosition" ] = {}      # initial translation
                interstellar[_id][ "initialQuaternion" ] = {}      # initial quaternion
                interstellar[_id][ "initialPosition" ][ "positionX" ] = translation[_id][ 0 ]      # initial translation
                interstellar[_id][ "initialPosition" ][ "positionY" ] = translation[_id][ 1 ]      # initial translation
                interstellar[_id][ "initialPosition" ][ "positionZ" ] = translation[_id][ 2 ]      # initial translation
                interstellar[_id][ "initialQuaternion" ][ "quaternionX" ] = orientation[_id][ 0 ]      # initial quaternion
                interstellar[_id][ "initialQuaternion" ][ "quaternionY" ] = orientation[_id][ 1 ]      # initial quaternion
                interstellar[_id][ "initialQuaternion" ][ "quaternionZ" ] = orientation[_id][ 2 ]      # initial quaternion
                interstellar[_id][ "initialQuaternion" ][ "quaternionW" ] = orientation[_id][ 3 ]      # initial quaternion
            if (_id < 0):
                interstellar[_id][ "positionX" ] = translation[_id][ 0 ]
                interstellar[_id][ "positionY" ] = translation[_id][ 1 ]
                interstellar[_id][ "positionZ" ] = translation[_id][ 2 ]
                interstellar[_id][ "quaternionX" ] = orientation[_id][ 1 ]
                interstellar[_id][ "quaternionY" ] = orientation[_id][ 2 ]
                interstellar[_id][ "quaternionZ" ] = orientation[_id][ 3 ]
                interstellar[_id][ "quaternionW" ] = orientation[_id][ 0 ]

        #print("interstellar:  %s" % interstellar )
        #print("interstellarFix:  %s" % interstellar )

        zeroIN = False                                          # we check if _id = 0 exists or not
        for instance in io.instances():
            _id = io.instances()[instance].attrs['id']
            if ( _id==0 ):
                zeroIN = True
                break


        for timestep in range(numberOfTimeSteps):
            #print("\nNumber of Time Steps: %d" % numberOfTimeSteps)
            update_progress(timestep/float(numberOfTimeSteps))
            positions = dpos_data[nbobjs*timestep:nbobjs*timestep+nbobjs, 2:]
            #print("\nPositions: %s" % positions)             #Positions of all objects at timestep
            for instance in io.instances():
                 _id = io.instances()[instance].attrs['id']
                 _idTab = _id
                 if ( _id >= 0 ):
                    if ( zeroIN ):
                        _idTab += 1                                                                 # if _id=0 present then we have to be careful that positions[0] returns positions of object id=0 whereas if _id=0 NOT PRESENT then position[0] returns positions of object id=1
                    q0, q1, q2, q3, q4, q5, q6 = [float(x) for x in positions[_idTab-1,:]]           # positions and quaternions of the object at timestep
                    #print("orientation["+str(_id)+"]",orientation[_id])
                                                                                                  # position coordinates and quaternions are gathered in interstellar
                    interstellar[_id][ "positionX" ].append(q0)
                    interstellar[_id][ "positionY" ].append(q1)
                    interstellar[_id][ "positionZ" ].append(q2)


                    interstellar[_id][ "quaternionW" ].append(q3)
                    interstellar[_id][ "quaternionX" ].append(q4)
                    interstellar[_id][ "quaternionY" ].append(q5)
                    interstellar[_id][ "quaternionZ" ].append(q6)
        my_ren = ThreejsRenderer(longueur=longueur, numberOfVectors = numberOfVectors, interstellar=interstellar, numberOfTimeSteps=numberOfTimeSteps, numberStaticObjects = numberStaticObjects, interstellarVectors = interstellarVectors )     #create the class with parameters
        my_ren.CreateDictionaryOfShapes(dictionaryOfShapes)                              # calling the method CreateDictionaryOfShapes to create the js file corresponding to the shapes


        def inplace_change(filename, old_string, new_string):
            s=open(filename).read()
            if old_string in s:
                print('Changing "{old_string}" to "{new_string}"'.format(**locals()), "in", filename)
                s=s.replace(old_string, new_string)
                f=open(filename, 'w')
                f.write(s)
                f.flush()
                f.close()
            else:
                print('No occurences of "{old_string}" found.'.format(**locals()))
        def insert_beginning_file(filename,  new_string):
            s=open(filename).read()
            print(s[0])
            s = new_string +s
            print(s[0])
            f=open(filename, 'w')
            f.write(s)
            f.flush()
            f.close()
           
            
        def insert_end_file(filename,  new_string):
            s=open(filename).read()
            print(s[-1])
            s = s +new_string 
            print(s[-1])
            f=open(filename, 'w')
            f.write(s)
            f.flush()
            f.close()
           
            
        #-----------Replace all "Shape" by "Shape1","Shape2",....etc...------------#
        for instance in io.instances():
            _id = io.instances()[instance].attrs['id']
            if ( _id >= 0 ):
                inplace_change(my_ren._path+"/shape"+str(_id)+".js", "Shape","Shape"+str(_id) ) # replace Shape in the Shape1 file by Shape1, same for Shape2, Shape3,...
            if ( _id < 0 ):
                inplace_change(my_ren._path+"/shape_"+str(abs(_id))+".js", "Shape","Shape_"+str(abs(_id)) ) # replace Shape in the Shape1 file by Shape1, same for Shape2, Shape3,...
             
        #---------------------Converting interstellar into JSON--------------------#
        import json
        writting = json.dumps(interstellar,ensure_ascii=False)                                                                                  # add "var dataInterstellar =" at the begining of the interstellar.json file and a "; at the end"
        with open(my_ren._path+'/interstellar.json', 'w') as outfile:
            json.dump(interstellar, outfile, indent=4)
        insert_beginning_file(os.path.join(my_ren._path, 'interstellar.json'),"var dataInterstellar =")
        insert_end_file(os.path.join(my_ren._path, 'interstellar.json'),";")


        writting = json.dumps(interstellarVectors,ensure_ascii=True)                                                                                  # add "var dataInterstellarVectors =" at the begining of the interstellarVectors.json file and a "; at the end"
        with open(my_ren._path+'/interstellarVectors.json', 'w') as outfile:
            json.dump(interstellarVectors, outfile, indent=4)
        insert_beginning_file(os.path.join(my_ren._path, 'interstellarVectors.json'),"var dataInterstellarVectors =")
        insert_end_file(os.path.join(my_ren._path, 'interstellarVectors.json'),";")

        #---Move the files to directory SliderCrank/Shape -------------------------#
        #subprocess.call(["mv " + my_ren._path+" "+os.getcwd() + "/simulation/static/renderer/Interstellar/shape_temps/"],shell=True)                 # move all the files from the TPF folder to the path indicated in green

"""
        urlstext = open(os.getcwd()+'/simulation/static/renderer/Interstellar/urls.txt', 'w')
        urlstext.write( my_ren._path[5:] )
        urlstext.close()
"""

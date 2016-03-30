#include "fc3d_NaturalMapABGenerated.h"
#include "funcodegen.h"
/*@
requires (0.0 <= rn <= 1.0e+6);
requires (-1.0e+6 <= rt1 <= 1.0e+6);
requires (-1.0e+6 <= rt2 <= 1.0e+6);
requires (-1.0e+6 <= un <= 1.0e+6);
requires (-1.0e+6 <= ut1 <= 1.0e+6);
requires (-1.0e+6 <= ut2 <= 1.0e+6);
requires (0.0 <= mu <= 1.0);
requires (-1.0e+6 <= rhon <= 1.0e+6);
requires (-1.0e+6 <= rhot1 <= 1.0e+6);
requires (-1.0e+6 <= rhot2 <= 1.0e+6);
assigns result[0..17];
ensures \is_finite((double) result[0]);
ensures \is_finite((double) result[1]);
ensures \is_finite((double) result[2]);
ensures \is_finite((double) result[3]);
ensures \is_finite((double) result[4]);
ensures \is_finite((double) result[5]);
ensures \is_finite((double) result[6]);
ensures \is_finite((double) result[7]);
ensures \is_finite((double) result[8]);
ensures \is_finite((double) result[9]);
ensures \is_finite((double) result[10]);
ensures \is_finite((double) result[11]);
ensures \is_finite((double) result[12]);
ensures \is_finite((double) result[13]);
ensures \is_finite((double) result[14]);
ensures \is_finite((double) result[15]);
ensures \is_finite((double) result[16]);
ensures \is_finite((double) result[17]);*/
void fc3d_NaturalMapABGenerated(
    double rn,
    double rt1,
    double rt2,
    double un,
    double ut1,
    double ut2,
    double mu,
    double rhon,
    double rhot1,
    double rhot2,
    double *result)
{
    /*@ assert \is_finite((double) rn); */
    /*@ assert \is_finite((double) ut1); */
    /*@ assert \is_finite((double) mu); */
    /*@ assert \is_finite((double) un); */
    /*@ assert \is_finite((double) rt1); */
    /*@ assert \is_finite((double) ut2); */
    /*@ assert \is_finite((double) rt2); */
    double x7 = 0.;
    double x8 = 0.;
    double x9 = 0.;
    double x10 = 0.;
    int x11 = 0;
    double x12 = 0.;
    double x13 = 0.;
    double x14 = 0.;
    double x15 = 0.;
    int x16 = 0;
    int x17 = 0;

    double x1 = 0.;
    double x2 = 0.;
    double x3 = 0.;
    double x4 = 0.;
    double x5 = 0.;
    double x6 = 0.;
    double x64 = 0.;
    double x65 = 0.;
    double x86 = 0.;
    double x109 = 0.;
    double x146 = 0.;
    double x147 = 0.;
    double x352 = 0.;
    double x353 = 0.;
    double x354 = 0.;
    double x355 = 0.;
    double x529 = 0.;


    /*@ assert (\is_finite((double) (rt1))); */
    x7 = rt1*rt1;
    /*@ assert \is_finite((double) (x7)); */
    /*@ assert (x7) >= 0; */

    /*@ assert (\is_finite((double) (rt2))); */
    x8 = rt2*rt2;
    /*@ assert \is_finite((double) (x8)); */
    /*@ assert (x8) >= 0; */
    x9 = x7 + x8;
    /*@ assert \is_finite((double) (x9)); */

    /*@ assert (\is_finite((double) (x9))); */
    /*@ assert (x9 >= 0); */
    x10 = sqrt(x9);
    /*@ assert \is_finite((double) (x10)); */
    /*@ assert (x10) >= 0; */
    /*@ assert (x10) > 2.22044604925e-16 ==> x9 > 4.930380657631323783823303533017413935457540219431394e-32; */
    x11 = x10 <= 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x11 <==> (x10 <= 0.0000000000000002220446049250313080847263336181640625); */

    /*@ assert (\is_finite((double) (ut1))); */
    x12 = ut1*ut1;
    /*@ assert \is_finite((double) (x12)); */
    /*@ assert (x12) >= 0; */

    /*@ assert (\is_finite((double) (ut2))); */
    x13 = ut2*ut2;
    /*@ assert \is_finite((double) (x13)); */
    /*@ assert (x13) >= 0; */
    x14 = x12 + x13;
    /*@ assert \is_finite((double) (x14)); */

    /*@ assert (\is_finite((double) (x14))); */
    /*@ assert (x14 >= 0); */
    x15 = sqrt(x14);
    /*@ assert \is_finite((double) (x15)); */
    /*@ assert (x15) >= 0; */
    /*@ assert (x15) > 2.22044604925e-16 ==> x14 > 4.930380657631323783823303533017413935457540219431394e-32; */
    x16 = x15 <= 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x16 <==> (x15 <= 0.0000000000000002220446049250313080847263336181640625); */
    x17 = x11 && x16;
    /*@ assert x17 <==> (x11 && x16); */

    int x27 = 0;
    int x28 = 0;

    double x18 = 0.;
    double x19 = 0.;
    double x20 = 0.;
    double x21 = 0.;
    double x22 = 0.;
    double x23 = 0.;
    double x24 = 0.;
    double x25 = 0.;
    double x26 = 0.;
    double x66 = 0.;
    double x67 = 0.;
    double x68 = 0.;
    double x69 = 0.;
    double x70 = 0.;
    double x95 = 0.;
    double x96 = 0.;
    double x97 = 0.;
    double x110 = 0.;
    double x111 = 0.;
    double x112 = 0.;
    double x117 = 0.;
    double x118 = 0.;
    double x119 = 0.;
    double x120 = 0.;
    double x121 = 0.;
    double x133 = 0.;
    double x134 = 0.;
    double x135 = 0.;

    x27 = x15 > 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x27 <==> (x15 > 0.0000000000000002220446049250313080847263336181640625); */
    x28 = x11 && x27;
    /*@ assert x28 <==> (x11 && x27); */

    int x37 = 0;
    int x38 = 0;
    int x39 = 0;
    double x40 = 0.;
    double x41 = 0.;
    double x42 = 0.;
    double x43 = 0.;
    double x44 = 0.;
    double x45 = 0.;
    double x46 = 0.;
    double x47 = 0.;
    double x48 = 0.;
    double x49 = 0.;
    int x50 = 0;
    int x51 = 0;
    int x52 = 0;
    int x53 = 0;

    double x29 = 0.;
    double x30 = 0.;
    double x31 = 0.;
    double x32 = 0.;
    double x33 = 0.;
    double x34 = 0.;
    double x35 = 0.;
    double x36 = 0.;
    double x71 = 0.;
    double x78 = 0.;
    double x113 = 0.;
    double x114 = 0.;
    double x122 = 0.;
    double x123 = 0.;
    double x124 = 0.;
    double x125 = 0.;
    double x144 = 0.;
    double x145 = 0.;
    double x148 = 0.;
    double x149 = 0.;
    double x150 = 0.;
    double x151 = 0.;
    double x152 = 0.;
    double x153 = 0.;
    double x170 = 0.;
    double x174 = 0.;
    double x175 = 0.;
    double x176 = 0.;
    double x177 = 0.;
    double x184 = 0.;
    double x185 = 0.;
    double x186 = 0.;
    double x191 = 0.;
    double x200 = 0.;
    double x201 = 0.;
    double x207 = 0.;
    double x210 = 0.;
    double x211 = 0.;
    double x236 = 0.;
    double x238 = 0.;
    double x240 = 0.;
    double x254 = 0.;
    double x257 = 0.;
    double x258 = 0.;
    double x259 = 0.;
    double x260 = 0.;
    double x261 = 0.;
    double x262 = 0.;
    double x263 = 0.;
    double x264 = 0.;
    double x294 = 0.;
    double x295 = 0.;
    double x296 = 0.;
    double x297 = 0.;
    double x298 = 0.;
    double x299 = 0.;
    double x300 = 0.;
    double x301 = 0.;
    double x302 = 0.;
    double x303 = 0.;
    double x304 = 0.;
    double x305 = 0.;
    double x306 = 0.;
    double x307 = 0.;
    double x308 = 0.;
    double x329 = 0.;
    double x330 = 0.;
    double x331 = 0.;
    double x340 = 0.;
    double x341 = 0.;
    double x342 = 0.;
    double x343 = 0.;
    double x344 = 0.;
    double x345 = 0.;
    double x346 = 0.;
    double x347 = 0.;
    double x348 = 0.;
    double x349 = 0.;
    double x350 = 0.;
    double x372 = 0.;
    double x377 = 0.;
    double x378 = 0.;
    double x380 = 0.;
    double x382 = 0.;
    double x383 = 0.;
    double x388 = 0.;
    double x390 = 0.;
    double x393 = 0.;
    double x395 = 0.;
    double x400 = 0.;
    double x414 = 0.;
    double x423 = 0.;
    double x428 = 0.;
    double x429 = 0.;
    double x430 = 0.;
    double x431 = 0.;
    double x433 = 0.;
    double x438 = 0.;
    double x440 = 0.;
    double x443 = 0.;
    double x458 = 0.;
    double x460 = 0.;
    double x462 = 0.;
    double x470 = 0.;
    double x482 = 0.;
    double x489 = 0.;
    double x500 = 0.;
    double x503 = 0.;
    double x521 = 0.;
    double x522 = 0.;
    double x523 = 0.;
    double x524 = 0.;
    double x525 = 0.;
    double x526 = 0.;
    double x530 = 0.;
    double x531 = 0.;
    double x532 = 0.;
    double x533 = 0.;
    double x534 = 0.;
    double x535 = 0.;
    double x536 = 0.;
    double x537 = 0.;
    double x538 = 0.;
    double x539 = 0.;
    double x540 = 0.;
    double x541 = 0.;
    double x542 = 0.;
    double x543 = 0.;
    double x544 = 0.;
    double x545 = 0.;
    double x546 = 0.;
    double x547 = 0.;
    double x548 = 0.;
    double x549 = 0.;
    double x550 = 0.;
    double x551 = 0.;
    double x562 = 0.;
    double x563 = 0.;
    double x564 = 0.;
    double x565 = 0.;
    double x566 = 0.;
    double x567 = 0.;
    double x568 = 0.;
    double x569 = 0.;
    double x570 = 0.;
    double x571 = 0.;
    double x572 = 0.;
    double x573 = 0.;
    double x574 = 0.;
    double x575 = 0.;
    double x576 = 0.;
    double x577 = 0.;
    double x578 = 0.;
    double x579 = 0.;
    double x580 = 0.;
    double x581 = 0.;
    double x582 = 0.;
    double x583 = 0.;
    double x584 = 0.;
    double x585 = 0.;

    x37 = x10 > 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x37 <==> (x10 > 0.0000000000000002220446049250313080847263336181640625); */
    x38 = x16 || x37;
    /*@ assert x38 <==> (x16 || x37); */
    x39 = x27 || x37;
    /*@ assert x39 <==> (x27 || x37); */
    x40 = -rt1;
    /*@ assert \is_finite((double) (x40)); */
    x41 = mu*ut1;
    /*@ assert \is_finite((double) (x41)); */
    x42 = x40 + x41;
    /*@ assert \is_finite((double) (x42)); */

    /*@ assert (\is_finite((double) (x42))); */
    x43 = x42*x42;
    /*@ assert \is_finite((double) (x43)); */
    /*@ assert (x43) >= 0; */
    x44 = -rt2;
    /*@ assert \is_finite((double) (x44)); */
    x45 = mu*ut2;
    /*@ assert \is_finite((double) (x45)); */
    x46 = x44 + x45;
    /*@ assert \is_finite((double) (x46)); */

    /*@ assert (\is_finite((double) (x46))); */
    x47 = x46*x46;
    /*@ assert \is_finite((double) (x47)); */
    /*@ assert (x47) >= 0; */
    x48 = x43 + x47;
    /*@ assert \is_finite((double) (x48)); */

    /*@ assert (\is_finite((double) (x48))); */
    /*@ assert (x48 >= 0); */
    x49 = sqrt(x48);
    /*@ assert \is_finite((double) (x49)); */
    /*@ assert (x49) >= 0; */
    /*@ assert (x49) > 2.22044604925e-16 ==> x48 > 4.930380657631323783823303533017413935457540219431394e-32; */
    x50 = x49 <= 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x50 <==> (x49 <= 0.0000000000000002220446049250313080847263336181640625); */
    x51 = x16 || x50;
    /*@ assert x51 <==> (x16 || x50); */
    x52 = x37 || x50;
    /*@ assert x52 <==> (x37 || x50); */
    x53 = x37 && x38 && x39 && x51 && x52;
    /*@ assert x53 <==> (x37 && x38 && x39 && x51 && x52); */

    int x63 = 0;

    double x54 = 0.;
    double x55 = 0.;
    double x56 = 0.;
    double x57 = 0.;
    double x58 = 0.;
    double x59 = 0.;
    double x60 = 0.;
    double x61 = 0.;
    double x62 = 0.;
    double x85 = 0.;
    double x87 = 0.;
    double x88 = 0.;
    double x89 = 0.;
    double x90 = 0.;
    double x91 = 0.;
    double x92 = 0.;
    double x93 = 0.;
    double x94 = 0.;
    double x102 = 0.;
    double x103 = 0.;
    double x104 = 0.;
    double x105 = 0.;
    double x106 = 0.;
    double x107 = 0.;
    double x108 = 0.;
    double x115 = 0.;
    double x116 = 0.;
    double x126 = 0.;
    double x127 = 0.;
    double x128 = 0.;
    double x129 = 0.;
    double x130 = 0.;
    double x131 = 0.;
    double x132 = 0.;
    double x136 = 0.;
    double x137 = 0.;
    double x138 = 0.;
    double x139 = 0.;
    double x140 = 0.;
    double x154 = 0.;
    double x213 = 0.;
    double x214 = 0.;
    double x215 = 0.;
    double x216 = 0.;
    double x217 = 0.;
    double x218 = 0.;
    double x219 = 0.;
    double x220 = 0.;
    double x221 = 0.;
    double x222 = 0.;
    double x223 = 0.;
    double x224 = 0.;
    double x225 = 0.;
    double x226 = 0.;
    double x227 = 0.;
    double x228 = 0.;
    double x266 = 0.;
    double x267 = 0.;
    double x268 = 0.;
    double x269 = 0.;
    double x270 = 0.;
    double x271 = 0.;
    double x272 = 0.;
    double x273 = 0.;
    double x276 = 0.;
    double x277 = 0.;
    double x278 = 0.;
    double x279 = 0.;
    double x309 = 0.;
    double x310 = 0.;
    double x311 = 0.;
    double x312 = 0.;
    double x313 = 0.;
    double x314 = 0.;
    double x315 = 0.;
    double x316 = 0.;
    double x317 = 0.;
    double x318 = 0.;
    double x319 = 0.;
    double x320 = 0.;
    double x321 = 0.;
    double x322 = 0.;
    double x323 = 0.;
    double x332 = 0.;
    double x333 = 0.;
    double x334 = 0.;
    double x335 = 0.;
    double x336 = 0.;
    double x337 = 0.;
    double x351 = 0.;
    double x447 = 0.;
    double x448 = 0.;
    double x449 = 0.;
    double x450 = 0.;
    double x511 = 0.;
    double x512 = 0.;
    double x513 = 0.;
    double x514 = 0.;
    double x515 = 0.;
    double x516 = 0.;
    double x517 = 0.;
    double x518 = 0.;
    double x527 = 0.;
    double x528 = 0.;
    double x552 = 0.;
    double x553 = 0.;
    double x586 = 0.;
    double x587 = 0.;
    double x588 = 0.;
    double x589 = 0.;
    double x590 = 0.;
    double x591 = 0.;
    double x592 = 0.;
    double x593 = 0.;
    double x594 = 0.;
    double x595 = 0.;

    x63 = x27 && x37 && x49 > 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x63 <==> (x27 && x37 && x49 > 0.0000000000000002220446049250313080847263336181640625); */

    int x82 = 0;

    double x72 = 0.;
    double x73 = 0.;
    double x74 = 0.;
    double x75 = 0.;
    double x76 = 0.;
    double x77 = 0.;
    double x79 = 0.;
    double x80 = 0.;
    double x81 = 0.;
    double x98 = 0.;
    double x99 = 0.;
    double x100 = 0.;
    double x101 = 0.;
    double x171 = 0.;
    double x172 = 0.;
    double x173 = 0.;
    double x178 = 0.;
    double x179 = 0.;
    double x180 = 0.;
    double x181 = 0.;
    double x182 = 0.;
    double x183 = 0.;
    double x187 = 0.;
    double x188 = 0.;
    double x189 = 0.;
    double x190 = 0.;
    double x192 = 0.;
    double x193 = 0.;
    double x194 = 0.;
    double x195 = 0.;
    double x196 = 0.;
    double x197 = 0.;
    double x198 = 0.;
    double x199 = 0.;
    double x202 = 0.;
    double x203 = 0.;
    double x204 = 0.;
    double x205 = 0.;
    double x235 = 0.;
    double x237 = 0.;
    double x239 = 0.;
    double x241 = 0.;
    double x242 = 0.;
    double x243 = 0.;
    double x244 = 0.;
    double x245 = 0.;
    double x246 = 0.;
    double x247 = 0.;
    double x248 = 0.;
    double x249 = 0.;
    double x250 = 0.;
    double x251 = 0.;
    double x252 = 0.;
    double x253 = 0.;
    double x255 = 0.;
    double x256 = 0.;
    double x356 = 0.;
    double x357 = 0.;
    double x358 = 0.;
    double x359 = 0.;
    double x360 = 0.;
    double x361 = 0.;
    double x362 = 0.;
    double x363 = 0.;
    double x364 = 0.;
    double x365 = 0.;
    double x366 = 0.;
    double x367 = 0.;
    double x368 = 0.;
    double x369 = 0.;
    double x370 = 0.;
    double x371 = 0.;
    double x373 = 0.;
    double x374 = 0.;
    double x375 = 0.;
    double x376 = 0.;
    double x379 = 0.;
    double x381 = 0.;
    double x384 = 0.;
    double x385 = 0.;
    double x386 = 0.;
    double x387 = 0.;
    double x389 = 0.;
    double x391 = 0.;
    double x392 = 0.;
    double x394 = 0.;
    double x396 = 0.;
    double x397 = 0.;
    double x398 = 0.;
    double x399 = 0.;
    double x401 = 0.;
    double x402 = 0.;
    double x403 = 0.;
    double x404 = 0.;
    double x405 = 0.;
    double x406 = 0.;
    double x407 = 0.;
    double x408 = 0.;
    double x409 = 0.;
    double x410 = 0.;
    double x411 = 0.;
    double x412 = 0.;
    double x413 = 0.;
    double x415 = 0.;
    double x416 = 0.;
    double x417 = 0.;
    double x418 = 0.;
    double x419 = 0.;
    double x420 = 0.;
    double x421 = 0.;
    double x422 = 0.;
    double x424 = 0.;
    double x425 = 0.;
    double x426 = 0.;
    double x459 = 0.;
    double x461 = 0.;
    double x463 = 0.;
    double x464 = 0.;
    double x465 = 0.;
    double x466 = 0.;
    double x467 = 0.;
    double x468 = 0.;
    double x469 = 0.;
    double x471 = 0.;
    double x472 = 0.;
    double x473 = 0.;
    double x474 = 0.;
    double x475 = 0.;
    double x476 = 0.;
    double x477 = 0.;
    double x478 = 0.;
    double x479 = 0.;
    double x480 = 0.;
    double x481 = 0.;
    double x483 = 0.;
    double x484 = 0.;
    double x485 = 0.;
    double x486 = 0.;
    double x487 = 0.;
    double x488 = 0.;
    double x490 = 0.;
    double x491 = 0.;

    x82 = x16 && x37;
    /*@ assert x82 <==> (x16 && x37); */

    int x84 = 0;

    double x83 = 0.;
    double x206 = 0.;
    double x208 = 0.;
    double x209 = 0.;
    double x212 = 0.;
    double x265 = 0.;
    double x427 = 0.;
    double x432 = 0.;
    double x434 = 0.;
    double x435 = 0.;
    double x436 = 0.;
    double x437 = 0.;
    double x439 = 0.;
    double x441 = 0.;
    double x442 = 0.;
    double x444 = 0.;
    double x445 = 0.;
    double x446 = 0.;
    double x492 = 0.;
    double x493 = 0.;
    double x494 = 0.;
    double x495 = 0.;
    double x496 = 0.;
    double x497 = 0.;
    double x498 = 0.;
    double x499 = 0.;
    double x501 = 0.;
    double x502 = 0.;
    double x504 = 0.;
    double x505 = 0.;
    double x506 = 0.;
    double x507 = 0.;
    double x508 = 0.;
    double x509 = 0.;
    double x510 = 0.;

    x84 = x27 && x37 && x50;
    /*@ assert x84 <==> (x27 && x37 && x50); */

    int x142 = 0;

    double x141 = 0.;
    double x155 = 0.;
    double x156 = 0.;
    double x157 = 0.;
    double x158 = 0.;
    double x159 = 0.;
    double x160 = 0.;
    double x161 = 0.;
    double x162 = 0.;
    double x163 = 0.;
    double x164 = 0.;
    double x165 = 0.;
    double x166 = 0.;
    double x167 = 0.;
    double x168 = 0.;
    double x169 = 0.;
    double x229 = 0.;
    double x230 = 0.;
    double x231 = 0.;
    double x232 = 0.;
    double x233 = 0.;
    double x234 = 0.;
    double x274 = 0.;
    double x275 = 0.;
    double x280 = 0.;
    double x281 = 0.;
    double x282 = 0.;
    double x283 = 0.;
    double x284 = 0.;
    double x285 = 0.;
    double x286 = 0.;
    double x287 = 0.;
    double x288 = 0.;
    double x289 = 0.;
    double x290 = 0.;
    double x291 = 0.;
    double x292 = 0.;
    double x293 = 0.;
    double x324 = 0.;
    double x325 = 0.;
    double x326 = 0.;
    double x327 = 0.;
    double x328 = 0.;
    double x338 = 0.;
    double x451 = 0.;
    double x452 = 0.;
    double x453 = 0.;
    double x454 = 0.;
    double x455 = 0.;
    double x456 = 0.;
    double x457 = 0.;
    double x519 = 0.;
    double x520 = 0.;
    double x554 = 0.;
    double x555 = 0.;
    double x556 = 0.;
    double x557 = 0.;
    double x558 = 0.;
    double x559 = 0.;
    double x560 = 0.;
    double x561 = 0.;

    x22 = mu*x15;
    /*@ assert \is_finite((double) (x22)); */
    x142 = x11 && x27 && x22 > 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x142 <==> (x11 && x27 && x22 > 0.0000000000000002220446049250313080847263336181640625); */

    int x143 = 0;

    double x339 = 0.;

    x143 = x11 && x27 && x22 <= 0.0000000000000002220446049250313080847263336181640625;
    /*@ assert x143 <==> (x11 && x27 && x22 <= 0.0000000000000002220446049250313080847263336181640625); */

    if (x17)
    {
        x1 = 1.0*un;
        /*@ assert \is_finite((double) (x1)); */
        x2 = -x1;
        /*@ assert \is_finite((double) (x2)); */
        x3 = mu*rn;
        /*@ assert \is_finite((double) (x3)); */
        x4 = 1.0*x3;
        /*@ assert \is_finite((double) (x4)); */
        x5 = x2 + x4;
        /*@ assert \is_finite((double) (x5)); */
        x6 = Heaviside(x5);
        /*@ assert \is_finite((double) (x6)); */
        x64 = 0.7071067811865475727373109293694142252206802368164063*mu;
        /*@ assert \is_finite((double) (x64)); */
        x65 = x6*x64;
        /*@ assert \is_finite((double) (x65)); */
        x86 = 1.0*mu;
        /*@ assert \is_finite((double) (x86)); */
        x109 = -x6*x86;
        /*@ assert \is_finite((double) (x109)); */
        x146 = -x4;
        /*@ assert \is_finite((double) (x146)); */
        x147 = x1 + x146;
        /*@ assert \is_finite((double) (x147)); */
        x352 = -x64;
        /*@ assert \is_finite((double) (x352)); */
        x353 = Heaviside(x147);
        /*@ assert \is_finite((double) (x353)); */
        x354 = x353*x64;
        /*@ assert \is_finite((double) (x354)); */
        x355 = x352 + x354;
        /*@ assert \is_finite((double) (x355)); */
        x529 = -0.7071067811865475727373109293694142252206802368164063*x353;
        /*@ assert \is_finite((double) (x529)); */

    }
    if (x28)
    {
        x3 = mu*rn;
        /*@ assert \is_finite((double) (x3)); */
        x18 = -un;
        /*@ assert \is_finite((double) (x18)); */
        x19 = x18 + x3;
        /*@ assert \is_finite((double) (x19)); */
        x20 = Heaviside(x19);
        /*@ assert \is_finite((double) (x20)); */
        x21 = 0.5*x20;
        /*@ assert \is_finite((double) (x21)); */
        x23 = -2*x22;
        /*@ assert \is_finite((double) (x23)); */
        x24 = x19 + x23;
        /*@ assert \is_finite((double) (x24)); */
        x25 = Heaviside(x24);
        /*@ assert \is_finite((double) (x25)); */
        x26 = 0.5*x25;
        /*@ assert \is_finite((double) (x26)); */
        x66 = mu*ut1;
        /*@ assert \is_finite((double) (x66)); */

        /*@ assert (\is_finite((double) (x15))); */
        /*@ assert (x15 < -1.09476442525e-47 || x15 > 1.09476442525e-47); */
        x67 = 1.0/x15;
        /*@ assert \is_finite((double) (x67)); */
        x68 = 1.0*x67;
        /*@ assert \is_finite((double) (x68)); */
        x69 = x66*x68;
        /*@ assert \is_finite((double) (x69)); */
        x70 = x25*x69;
        /*@ assert \is_finite((double) (x70)); */
        x95 = mu*ut2;
        /*@ assert \is_finite((double) (x95)); */
        x96 = x68*x95;
        /*@ assert \is_finite((double) (x96)); */
        x97 = x25*x96;
        /*@ assert \is_finite((double) (x97)); */
        x110 = -mu*x21;
        /*@ assert \is_finite((double) (x110)); */
        x111 = mu*x26;
        /*@ assert \is_finite((double) (x111)); */
        x112 = -x111;
        /*@ assert \is_finite((double) (x112)); */
        x117 = 0.5*x20*x67;
        /*@ assert \is_finite((double) (x117)); */
        x118 = ut1*x117;
        /*@ assert \is_finite((double) (x118)); */
        x119 = 0.5*x25*x67;
        /*@ assert \is_finite((double) (x119)); */
        x120 = ut1*x119;
        /*@ assert \is_finite((double) (x120)); */
        x121 = -x120;
        /*@ assert \is_finite((double) (x121)); */
        x133 = ut2*x117;
        /*@ assert \is_finite((double) (x133)); */
        x134 = ut2*x119;
        /*@ assert \is_finite((double) (x134)); */
        x135 = -x134;
        /*@ assert \is_finite((double) (x135)); */

    }
    if (x53)
    {
        x1 = 1.0*un;
        /*@ assert \is_finite((double) (x1)); */
        x2 = -x1;
        /*@ assert \is_finite((double) (x2)); */
        x3 = mu*rn;
        /*@ assert \is_finite((double) (x3)); */
        x4 = 1.0*x3;
        /*@ assert \is_finite((double) (x4)); */
        x5 = x2 + x4;
        /*@ assert \is_finite((double) (x5)); */
        x18 = -un;
        /*@ assert \is_finite((double) (x18)); */
        x29 = 1.0*x10;
        /*@ assert \is_finite((double) (x29)); */
        x30 = -x29;
        /*@ assert \is_finite((double) (x30)); */
        x31 = x30 + x5;
        /*@ assert \is_finite((double) (x31)); */
        x32 = Heaviside(x31);
        /*@ assert \is_finite((double) (x32)); */
        x33 = 0.5*x32;
        /*@ assert \is_finite((double) (x33)); */
        x34 = x29 + x5;
        /*@ assert \is_finite((double) (x34)); */
        x35 = Heaviside(x34);
        /*@ assert \is_finite((double) (x35)); */
        x36 = 0.5*x35;
        /*@ assert \is_finite((double) (x36)); */

        /*@ assert (\is_finite((double) (x10))); */
        /*@ assert (x10 < -1.09476442525e-47 || x10 > 1.09476442525e-47); */
        x71 = 1.0/x10;
        /*@ assert \is_finite((double) (x71)); */
        x78 = x10*x35;
        /*@ assert \is_finite((double) (x78)); */
        x113 = -mu*x33;
        /*@ assert \is_finite((double) (x113)); */
        x114 = -mu*x36;
        /*@ assert \is_finite((double) (x114)); */
        x122 = 1.0*x35;
        /*@ assert \is_finite((double) (x122)); */
        x123 = -x122;
        /*@ assert \is_finite((double) (x123)); */
        x124 = x123 + x32;
        /*@ assert \is_finite((double) (x124)); */
        x125 = 0.5*x124*x71;
        /*@ assert \is_finite((double) (x125)); */
        x144 = Heaviside(x10);
        /*@ assert \is_finite((double) (x144)); */
        x145 = 0.5*x144;
        /*@ assert \is_finite((double) (x145)); */
        x146 = -x4;
        /*@ assert \is_finite((double) (x146)); */
        x147 = x1 + x146;
        /*@ assert \is_finite((double) (x147)); */
        x148 = Heaviside(x147 + x30);
        /*@ assert \is_finite((double) (x148)); */
        x149 = Heaviside(x147 + x29);
        /*@ assert \is_finite((double) (x149)); */
        x150 = 1.0*x149;
        /*@ assert \is_finite((double) (x150)); */
        x151 = -x150;
        /*@ assert \is_finite((double) (x151)); */
        x152 = x148 + x151;
        /*@ assert \is_finite((double) (x152)); */
        x153 = rt1*x145*x152*x71;
        /*@ assert \is_finite((double) (x153)); */

        /*@ assert (\is_finite((double) (x9))); */
        /*@ assert (x9 < -1.09476442525e-47 || x9 > 1.09476442525e-47); */
        x170 = 1.0/(x9*x9);
        /*@ assert \is_finite((double) (x170)); */
        /*@ assert (x170) >= 0; */

        /*@ assert (\is_finite((double) (rt2))); */
        x174 = rt2*rt2*rt2*rt2;
        /*@ assert \is_finite((double) (x174)); */
        x175 = mu*rn;
        /*@ assert \is_finite((double) (x175)); */
        x176 = x175 + x18;
        /*@ assert \is_finite((double) (x176)); */
        x177 = Max(0, x10 + x176);
        /*@ assert \is_finite((double) (x177)); */
        /*@ assert (x177) >= 0; */
        x184 = -x10;
        /*@ assert \is_finite((double) (x184)); */
        x185 = Max(0, x176 + x184);
        /*@ assert \is_finite((double) (x185)); */
        /*@ assert (x185) >= 0; */
        x186 = x174*x185;
        /*@ assert \is_finite((double) (x186)); */
        x191 = x177*x7*x8;
        /*@ assert \is_finite((double) (x191)); */
        x200 = x10*x9;
        /*@ assert \is_finite((double) (x200)); */
        x201 = 2.0*x200;
        /*@ assert \is_finite((double) (x201)); */
        x207 = x174*x177;
        /*@ assert \is_finite((double) (x207)); */
        x210 = x185*x7;
        /*@ assert \is_finite((double) (x210)); */
        x211 = x210*x8;
        /*@ assert \is_finite((double) (x211)); */

        /*@ assert (\is_finite((double) (rt2))); */
        x236 = rt2*rt2*rt2;
        /*@ assert \is_finite((double) (x236)); */
        x238 = x177*x7;
        /*@ assert \is_finite((double) (x238)); */

        /*@ assert (\is_finite((double) (rt1))); */
        x240 = rt1*rt1*rt1*rt1;
        /*@ assert \is_finite((double) (x240)); */
        x254 = rt2*x149;
        /*@ assert \is_finite((double) (x254)); */
        x257 = 0.5*rt1*rt2*x144;
        /*@ assert \is_finite((double) (x257)); */
        x258 = -x201;
        /*@ assert \is_finite((double) (x258)); */
        x259 = x177*x8;
        /*@ assert \is_finite((double) (x259)); */
        x260 = 1.0*x185;
        /*@ assert \is_finite((double) (x260)); */
        x261 = x260*x7;
        /*@ assert \is_finite((double) (x261)); */
        x262 = -x261;
        /*@ assert \is_finite((double) (x262)); */
        x263 = -x260*x8;
        /*@ assert \is_finite((double) (x263)); */
        x264 = x238 + x258 + x259 + x262 + x263;
        /*@ assert \is_finite((double) (x264)); */
        x294 = 1.0*x170*x71;
        /*@ assert \is_finite((double) (x294)); */
        x295 = x200*x7;
        /*@ assert \is_finite((double) (x295)); */
        x296 = x200*x8;
        /*@ assert \is_finite((double) (x296)); */
        x297 = x295 + x296;
        /*@ assert \is_finite((double) (x297)); */
        x298 = 1.0*x10*x144*x9;
        /*@ assert \is_finite((double) (x298)); */
        x299 = -x298*x7;
        /*@ assert \is_finite((double) (x299)); */
        x300 = -x145*x207;
        /*@ assert \is_finite((double) (x300)); */
        x301 = x145*x186;
        /*@ assert \is_finite((double) (x301)); */
        x302 = -x145*x191;
        /*@ assert \is_finite((double) (x302)); */
        x303 = x145*x211;
        /*@ assert \is_finite((double) (x303)); */
        x304 = x10*x148*x7*x9;
        /*@ assert \is_finite((double) (x304)); */
        x305 = x145*x304;
        /*@ assert \is_finite((double) (x305)); */
        x306 = 0.5*x144*x149;
        /*@ assert \is_finite((double) (x306)); */
        x307 = x295*x306;
        /*@ assert \is_finite((double) (x307)); */
        x308 = x297 + x299 + x300 + x301 + x302 + x303 + x305 + x307;
        /*@ assert \is_finite((double) (x308)); */
        x329 = x148*x200;
        /*@ assert \is_finite((double) (x329)); */
        x330 = x149*x200;
        /*@ assert \is_finite((double) (x330)); */
        x331 = x264 + x329 + x330;
        /*@ assert \is_finite((double) (x331)); */
        x340 = -x10*x33;
        /*@ assert \is_finite((double) (x340)); */
        x341 = x10*x36;
        /*@ assert \is_finite((double) (x341)); */
        x342 = rt2*x144*x148;
        /*@ assert \is_finite((double) (x342)); */
        x343 = -0.5*x342;
        /*@ assert \is_finite((double) (x343)); */
        x344 = rt2*x144;
        /*@ assert \is_finite((double) (x344)); */
        x345 = x149*x344;
        /*@ assert \is_finite((double) (x345)); */
        x346 = 0.5*x345;
        /*@ assert \is_finite((double) (x346)); */
        x347 = x10*x148;
        /*@ assert \is_finite((double) (x347)); */
        x348 = x145*x347;
        /*@ assert \is_finite((double) (x348)); */
        x349 = -x10*x306;
        /*@ assert \is_finite((double) (x349)); */
        x350 = x340 + x341 + x343 + x346 + x348 + x349;
        /*@ assert \is_finite((double) (x350)); */
        x372 = rt2*x10*x9;
        /*@ assert \is_finite((double) (x372)); */

        /*@ assert (\is_finite((double) (rt2))); */
        x377 = rt2*rt2*rt2*rt2*rt2;
        /*@ assert \is_finite((double) (x377)); */
        x378 = x144*x148*x377;
        /*@ assert \is_finite((double) (x378)); */
        x380 = x144*x149*x377;
        /*@ assert \is_finite((double) (x380)); */
        x382 = 2.0*x144;
        /*@ assert \is_finite((double) (x382)); */
        x383 = x174*x382;
        /*@ assert \is_finite((double) (x383)); */
        x388 = rt2*x144*x148*x240;
        /*@ assert \is_finite((double) (x388)); */
        x390 = rt2*x144*x149*x240;
        /*@ assert \is_finite((double) (x390)); */
        x393 = Max(0, x31);
        /*@ assert \is_finite((double) (x393)); */
        /*@ assert (x393) >= 0; */
        x395 = Max(0, x34);
        /*@ assert \is_finite((double) (x395)); */
        /*@ assert (x395) >= 0; */
        x400 = x144*x148*x236*x7;
        /*@ assert \is_finite((double) (x400)); */
        x414 = x372*x382;
        /*@ assert \is_finite((double) (x414)); */
        x423 = x10*x148*x8*x9;
        /*@ assert \is_finite((double) (x423)); */
        x428 = x240*x382;
        /*@ assert \is_finite((double) (x428)); */
        x429 = 4.0*x144*x7;
        /*@ assert \is_finite((double) (x429)); */
        x430 = x429*x8;
        /*@ assert \is_finite((double) (x430)); */
        x431 = x240*x35;
        /*@ assert \is_finite((double) (x431)); */
        x433 = x174*x35;
        /*@ assert \is_finite((double) (x433)); */
        x438 = x144*x236;
        /*@ assert \is_finite((double) (x438)); */
        x440 = x395*x438;
        /*@ assert \is_finite((double) (x440)); */
        x443 = rt2*x144*x7;
        /*@ assert \is_finite((double) (x443)); */
        x458 = x144*x377;
        /*@ assert \is_finite((double) (x458)); */
        x460 = rt2*x240;
        /*@ assert \is_finite((double) (x460)); */
        x462 = x236*x7;
        /*@ assert \is_finite((double) (x462)); */
        x470 = x236*x32*x7;
        /*@ assert \is_finite((double) (x470)); */
        x482 = 2.0*x144*x7*x8;
        /*@ assert \is_finite((double) (x482)); */
        x489 = x144*x149;
        /*@ assert \is_finite((double) (x489)); */
        x500 = x144*x240;
        /*@ assert \is_finite((double) (x500)); */
        x503 = 1.0*x395;
        /*@ assert \is_finite((double) (x503)); */
        x521 = 0.5*x71;
        /*@ assert \is_finite((double) (x521)); */
        x522 = -x29*x32;
        /*@ assert \is_finite((double) (x522)); */
        x523 = -1.0*x342;
        /*@ assert \is_finite((double) (x523)); */
        x524 = x144*x347;
        /*@ assert \is_finite((double) (x524)); */
        x525 = -x29*x489;
        /*@ assert \is_finite((double) (x525)); */
        x526 = x345 + x522 + x523 + x524 + x525 + x78;
        /*@ assert \is_finite((double) (x526)); */
        x530 = -x428;
        /*@ assert \is_finite((double) (x530)); */
        x531 = -x383;
        /*@ assert \is_finite((double) (x531)); */
        x532 = -x430;
        /*@ assert \is_finite((double) (x532)); */
        x533 = x240*x32;
        /*@ assert \is_finite((double) (x533)); */
        x534 = x174*x32;
        /*@ assert \is_finite((double) (x534)); */
        x535 = 2.0*x7*x8;
        /*@ assert \is_finite((double) (x535)); */
        x536 = x32*x535;
        /*@ assert \is_finite((double) (x536)); */
        x537 = x35*x535;
        /*@ assert \is_finite((double) (x537)); */
        x538 = x148*x500;
        /*@ assert \is_finite((double) (x538)); */
        x539 = x149*x500;
        /*@ assert \is_finite((double) (x539)); */
        x540 = x144*x174;
        /*@ assert \is_finite((double) (x540)); */
        x541 = x148*x540;
        /*@ assert \is_finite((double) (x541)); */
        x542 = x149*x540;
        /*@ assert \is_finite((double) (x542)); */
        x543 = x393*x438;
        /*@ assert \is_finite((double) (x543)); */
        x544 = -1.0*x440;
        /*@ assert \is_finite((double) (x544)); */
        x545 = x393*x443;
        /*@ assert \is_finite((double) (x545)); */
        x546 = -x443*x503;
        /*@ assert \is_finite((double) (x546)); */
        x547 = x148*x482;
        /*@ assert \is_finite((double) (x547)); */
        x548 = x149*x482;
        /*@ assert \is_finite((double) (x548)); */
        x549 = -rt2*x148*x298;
        /*@ assert \is_finite((double) (x549)); */
        x550 = -x254*x298;
        /*@ assert \is_finite((double) (x550)); */
        x551 = x414 + x431 + x433 + x530 + x531 + x532 + x533 + x534 + x536 + x537 + x538 + x539 + x541 + x542 + x543 + x544 + x545 + x546 + x547 + x548 + x549 + x550;
        /*@ assert \is_finite((double) (x551)); */
        x562 = x144*x460;
        /*@ assert \is_finite((double) (x562)); */
        x563 = x382*x462;
        /*@ assert \is_finite((double) (x563)); */
        x564 = -x33*x377;
        /*@ assert \is_finite((double) (x564)); */
        x565 = -x36*x377;
        /*@ assert \is_finite((double) (x565)); */
        x566 = -x33*x460;
        /*@ assert \is_finite((double) (x566)); */
        x567 = -x36*x460;
        /*@ assert \is_finite((double) (x567)); */
        x568 = -1.0*x470;
        /*@ assert \is_finite((double) (x568)); */
        x569 = -x122*x462;
        /*@ assert \is_finite((double) (x569)); */
        x570 = -x298*x8;
        /*@ assert \is_finite((double) (x570)); */
        x571 = -0.5*x378;
        /*@ assert \is_finite((double) (x571)); */
        x572 = -0.5*x380;
        /*@ assert \is_finite((double) (x572)); */
        x573 = 0.5*x144*x240;
        /*@ assert \is_finite((double) (x573)); */
        x574 = x393*x573;
        /*@ assert \is_finite((double) (x574)); */
        x575 = -x395*x573;
        /*@ assert \is_finite((double) (x575)); */
        x576 = -0.5*x388;
        /*@ assert \is_finite((double) (x576)); */
        x577 = -0.5*x390;
        /*@ assert \is_finite((double) (x577)); */
        x578 = -1.0*x400;
        /*@ assert \is_finite((double) (x578)); */
        x579 = -x144*x150*x236*x7;
        /*@ assert \is_finite((double) (x579)); */
        x580 = 0.5*x144*x7*x8;
        /*@ assert \is_finite((double) (x580)); */
        x581 = x393*x580;
        /*@ assert \is_finite((double) (x581)); */
        x582 = -x395*x580;
        /*@ assert \is_finite((double) (x582)); */
        x583 = x145*x423;
        /*@ assert \is_finite((double) (x583)); */
        x584 = x296*x306;
        /*@ assert \is_finite((double) (x584)); */
        x585 = x297 + x458 + x562 + x563 + x564 + x565 + x566 + x567 + x568 + x569 + x570 + x571 + x572 + x574 + x575 + x576 + x577 + x578 + x579 + x581 + x582 + x583 + x584;
        /*@ assert \is_finite((double) (x585)); */

    }
    if (x63)
    {
        x3 = mu*rn;
        /*@ assert \is_finite((double) (x3)); */
        x18 = -un;
        /*@ assert \is_finite((double) (x18)); */
        x54 = -x22;
        /*@ assert \is_finite((double) (x54)); */
        x55 = x18 + x3 + x54;
        /*@ assert \is_finite((double) (x55)); */
        x56 = x49 + x55;
        /*@ assert \is_finite((double) (x56)); */
        x57 = Heaviside(x56);
        /*@ assert \is_finite((double) (x57)); */
        x58 = 0.5*x57;
        /*@ assert \is_finite((double) (x58)); */
        x59 = -x49;
        /*@ assert \is_finite((double) (x59)); */
        x60 = x55 + x59;
        /*@ assert \is_finite((double) (x60)); */
        x61 = Heaviside(x60);
        /*@ assert \is_finite((double) (x61)); */
        x62 = 0.5*x61;
        /*@ assert \is_finite((double) (x62)); */
        x66 = mu*ut1;
        /*@ assert \is_finite((double) (x66)); */

        /*@ assert (\is_finite((double) (x15))); */
        /*@ assert (x15 < -1.09476442525e-47 || x15 > 1.09476442525e-47); */
        x67 = 1.0/x15;
        /*@ assert \is_finite((double) (x67)); */
        x68 = 1.0*x67;
        /*@ assert \is_finite((double) (x68)); */
        x69 = x66*x68;
        /*@ assert \is_finite((double) (x69)); */
        x85 = -x69;
        /*@ assert \is_finite((double) (x85)); */
        x86 = 1.0*mu;
        /*@ assert \is_finite((double) (x86)); */

        /*@ assert (\is_finite((double) (x49))); */
        /*@ assert (x49 < -1.09476442525e-47 || x49 > 1.09476442525e-47); */
        x87 = 1.0/x49;
        /*@ assert \is_finite((double) (x87)); */
        x88 = x86*x87;
        /*@ assert \is_finite((double) (x88)); */
        x89 = x42*x88;
        /*@ assert \is_finite((double) (x89)); */
        x90 = x85 + x89;
        /*@ assert \is_finite((double) (x90)); */
        x91 = -x58*x90;
        /*@ assert \is_finite((double) (x91)); */
        x92 = -x89;
        /*@ assert \is_finite((double) (x92)); */
        x93 = x85 + x92;
        /*@ assert \is_finite((double) (x93)); */
        x94 = -x62*x93;
        /*@ assert \is_finite((double) (x94)); */
        x95 = mu*ut2;
        /*@ assert \is_finite((double) (x95)); */
        x96 = x68*x95;
        /*@ assert \is_finite((double) (x96)); */
        x102 = -x96;
        /*@ assert \is_finite((double) (x102)); */
        x103 = x46*x88;
        /*@ assert \is_finite((double) (x103)); */
        x104 = x102 + x103;
        /*@ assert \is_finite((double) (x104)); */
        x105 = -x104*x58;
        /*@ assert \is_finite((double) (x105)); */
        x106 = -x103;
        /*@ assert \is_finite((double) (x106)); */
        x107 = x102 + x106;
        /*@ assert \is_finite((double) (x107)); */
        x108 = -x107*x62;
        /*@ assert \is_finite((double) (x108)); */
        x115 = -mu*x58;
        /*@ assert \is_finite((double) (x115)); */
        x116 = -mu*x62;
        /*@ assert \is_finite((double) (x116)); */
        x126 = -x66;
        /*@ assert \is_finite((double) (x126)); */
        x127 = rt1 + x126;
        /*@ assert \is_finite((double) (x127)); */
        x128 = 0.5*x57*x87;
        /*@ assert \is_finite((double) (x128)); */
        x129 = x127*x128;
        /*@ assert \is_finite((double) (x129)); */
        x130 = -x129;
        /*@ assert \is_finite((double) (x130)); */
        x131 = 0.5*x61*x87;
        /*@ assert \is_finite((double) (x131)); */
        x132 = x127*x131;
        /*@ assert \is_finite((double) (x132)); */
        x136 = -x95;
        /*@ assert \is_finite((double) (x136)); */
        x137 = rt2 + x136;
        /*@ assert \is_finite((double) (x137)); */
        x138 = x128*x137;
        /*@ assert \is_finite((double) (x138)); */
        x139 = -x138;
        /*@ assert \is_finite((double) (x139)); */
        x140 = x131*x137;
        /*@ assert \is_finite((double) (x140)); */
        x154 = x131*x42;
        /*@ assert \is_finite((double) (x154)); */
        x213 = -x129*x90;
        /*@ assert \is_finite((double) (x213)); */
        x214 = -x154*x93;
        /*@ assert \is_finite((double) (x214)); */
        x215 = Max(0, x60);
        /*@ assert \is_finite((double) (x215)); */
        /*@ assert (x215) >= 0; */
        x216 = 0.5*x215;
        /*@ assert \is_finite((double) (x216)); */

        /*@ assert (\is_finite((double) (x48))); */
        /*@ assert (x48 < -1.09476442525e-47 || x48 > 1.09476442525e-47); */
        x217 = 1.0/x48;
        /*@ assert \is_finite((double) (x217)); */
        x218 = 1.0*mu*x217*x87;
        /*@ assert \is_finite((double) (x218)); */
        x219 = -x218*x43;
        /*@ assert \is_finite((double) (x219)); */
        x220 = x219 + x88;
        /*@ assert \is_finite((double) (x220)); */
        x221 = -x216*x220;
        /*@ assert \is_finite((double) (x221)); */
        x222 = Max(0, x56);
        /*@ assert \is_finite((double) (x222)); */
        /*@ assert (x222) >= 0; */
        x223 = 0.5*x222;
        /*@ assert \is_finite((double) (x223)); */
        x224 = -x88;
        /*@ assert \is_finite((double) (x224)); */
        x225 = x127*x217;
        /*@ assert \is_finite((double) (x225)); */
        x226 = -x225*x89;
        /*@ assert \is_finite((double) (x226)); */
        x227 = x224 + x226;
        /*@ assert \is_finite((double) (x227)); */
        x228 = -x223*x227;
        /*@ assert \is_finite((double) (x228)); */
        x266 = x127*x217*x46;
        /*@ assert \is_finite((double) (x266)); */
        x267 = 0.5*x222*x87;
        /*@ assert \is_finite((double) (x267)); */
        x268 = x266*x267;
        /*@ assert \is_finite((double) (x268)); */
        x269 = mu*x268;
        /*@ assert \is_finite((double) (x269)); */
        x270 = 0.5*x215*x217*x42*x46*x87;
        /*@ assert \is_finite((double) (x270)); */
        x271 = mu*x270;
        /*@ assert \is_finite((double) (x271)); */
        x272 = -x104*x129;
        /*@ assert \is_finite((double) (x272)); */
        x273 = -x107*x154;
        /*@ assert \is_finite((double) (x273)); */
        x276 = 0.5*mu*x57*x87;
        /*@ assert \is_finite((double) (x276)); */
        x277 = -x127*x276;
        /*@ assert \is_finite((double) (x277)); */
        x278 = 0.5*mu*x61*x87;
        /*@ assert \is_finite((double) (x278)); */
        x279 = -x278*x42;
        /*@ assert \is_finite((double) (x279)); */

        /*@ assert (\is_finite((double) (x127))); */
        x309 = x127*x127;
        /*@ assert \is_finite((double) (x309)); */
        /*@ assert (x309) >= 0; */
        x310 = 0.5*x217*x57;
        /*@ assert \is_finite((double) (x310)); */
        x311 = -x309*x310;
        /*@ assert \is_finite((double) (x311)); */
        x312 = x127*x42;
        /*@ assert \is_finite((double) (x312)); */
        x313 = 0.5*x217*x61;
        /*@ assert \is_finite((double) (x313)); */
        x314 = x312*x313;
        /*@ assert \is_finite((double) (x314)); */
        x315 = 1.0*x87;
        /*@ assert \is_finite((double) (x315)); */
        x316 = -x315;
        /*@ assert \is_finite((double) (x316)); */
        x317 = 1.0*x217*x87;
        /*@ assert \is_finite((double) (x317)); */
        x318 = x317*x43;
        /*@ assert \is_finite((double) (x318)); */
        x319 = x316 + x318;
        /*@ assert \is_finite((double) (x319)); */
        x320 = -x216*x319;
        /*@ assert \is_finite((double) (x320)); */
        x321 = x312*x317;
        /*@ assert \is_finite((double) (x321)); */
        x322 = x315 + x321;
        /*@ assert \is_finite((double) (x322)); */
        x323 = -x223*x322;
        /*@ assert \is_finite((double) (x323)); */
        x332 = -0.5*x137*x225*x57;
        /*@ assert \is_finite((double) (x332)); */
        x333 = -x270;
        /*@ assert \is_finite((double) (x333)); */
        x334 = x332 + x333;
        /*@ assert \is_finite((double) (x334)); */
        x335 = x137*x217*x42;
        /*@ assert \is_finite((double) (x335)); */
        x336 = x335*x62;
        /*@ assert \is_finite((double) (x336)); */
        x337 = -x268;
        /*@ assert \is_finite((double) (x337)); */
        x351 = x131*x46;
        /*@ assert \is_finite((double) (x351)); */
        x447 = x267*x335;
        /*@ assert \is_finite((double) (x447)); */
        x448 = mu*x447;
        /*@ assert \is_finite((double) (x448)); */
        x449 = -x138*x90;
        /*@ assert \is_finite((double) (x449)); */
        x450 = -x351*x93;
        /*@ assert \is_finite((double) (x450)); */
        x511 = -x104*x138;
        /*@ assert \is_finite((double) (x511)); */
        x512 = -x107*x351;
        /*@ assert \is_finite((double) (x512)); */
        x513 = -x218*x47;
        /*@ assert \is_finite((double) (x513)); */
        x514 = x513 + x88;
        /*@ assert \is_finite((double) (x514)); */
        x515 = -x216*x514;
        /*@ assert \is_finite((double) (x515)); */
        x516 = -x103*x137*x217;
        /*@ assert \is_finite((double) (x516)); */
        x517 = x224 + x516;
        /*@ assert \is_finite((double) (x517)); */
        x518 = -x223*x517;
        /*@ assert \is_finite((double) (x518)); */
        x527 = -x137*x276;
        /*@ assert \is_finite((double) (x527)); */
        x528 = -x278*x46;
        /*@ assert \is_finite((double) (x528)); */
        x552 = x266*x62;
        /*@ assert \is_finite((double) (x552)); */
        x553 = -x447;
        /*@ assert \is_finite((double) (x553)); */

        /*@ assert (\is_finite((double) (x137))); */
        x586 = x137*x137;
        /*@ assert \is_finite((double) (x586)); */
        /*@ assert (x586) >= 0; */
        x587 = -x310*x586;
        /*@ assert \is_finite((double) (x587)); */
        x588 = x137*x46;
        /*@ assert \is_finite((double) (x588)); */
        x589 = x313*x588;
        /*@ assert \is_finite((double) (x589)); */
        x590 = x317*x47;
        /*@ assert \is_finite((double) (x590)); */
        x591 = x316 + x590;
        /*@ assert \is_finite((double) (x591)); */
        x592 = -x216*x591;
        /*@ assert \is_finite((double) (x592)); */
        x593 = x317*x588;
        /*@ assert \is_finite((double) (x593)); */
        x594 = x315 + x593;
        /*@ assert \is_finite((double) (x594)); */
        x595 = -x223*x594;
        /*@ assert \is_finite((double) (x595)); */

    }
    if (x82)
    {
        x1 = 1.0*un;
        /*@ assert \is_finite((double) (x1)); */
        x2 = -x1;
        /*@ assert \is_finite((double) (x2)); */
        x3 = mu*rn;
        /*@ assert \is_finite((double) (x3)); */
        x4 = 1.0*x3;
        /*@ assert \is_finite((double) (x4)); */
        x5 = x2 + x4;
        /*@ assert \is_finite((double) (x5)); */
        x18 = -un;
        /*@ assert \is_finite((double) (x18)); */
        x29 = 1.0*x10;
        /*@ assert \is_finite((double) (x29)); */
        x30 = -x29;
        /*@ assert \is_finite((double) (x30)); */
        x31 = x30 + x5;
        /*@ assert \is_finite((double) (x31)); */
        x32 = Heaviside(x31);
        /*@ assert \is_finite((double) (x32)); */
        x34 = x29 + x5;
        /*@ assert \is_finite((double) (x34)); */
        x35 = Heaviside(x34);
        /*@ assert \is_finite((double) (x35)); */

        /*@ assert (\is_finite((double) (x10))); */
        /*@ assert (x10 < -1.09476442525e-47 || x10 > 1.09476442525e-47); */
        x71 = 1.0/x10;
        /*@ assert \is_finite((double) (x71)); */
        x72 = 0.25*mu*x71;
        /*@ assert \is_finite((double) (x72)); */
        x73 = 2.0*rt1;
        /*@ assert \is_finite((double) (x73)); */
        x74 = -x32*x73;
        /*@ assert \is_finite((double) (x74)); */
        x75 = x35*x73;
        /*@ assert \is_finite((double) (x75)); */
        x76 = 1.414213562373095145474621858738828450441360473632813*x32;
        /*@ assert \is_finite((double) (x76)); */
        x77 = x10*x76;
        /*@ assert \is_finite((double) (x77)); */
        x78 = x10*x35;
        /*@ assert \is_finite((double) (x78)); */
        x79 = 1.414213562373095145474621858738828450441360473632813*x78;
        /*@ assert \is_finite((double) (x79)); */
        x80 = x77 + x79;
        /*@ assert \is_finite((double) (x80)); */
        x81 = x74 + x75 + x80;
        /*@ assert \is_finite((double) (x81)); */
        x98 = 2.0*rt2;
        /*@ assert \is_finite((double) (x98)); */
        x99 = -x32*x98;
        /*@ assert \is_finite((double) (x99)); */
        x100 = x35*x98;
        /*@ assert \is_finite((double) (x100)); */
        x101 = x100 + x80 + x99;
        /*@ assert \is_finite((double) (x101)); */
        x144 = Heaviside(x10);
        /*@ assert \is_finite((double) (x144)); */
        x146 = -x4;
        /*@ assert \is_finite((double) (x146)); */
        x147 = x1 + x146;
        /*@ assert \is_finite((double) (x147)); */
        x148 = Heaviside(x147 + x30);
        /*@ assert \is_finite((double) (x148)); */
        x149 = Heaviside(x147 + x29);
        /*@ assert \is_finite((double) (x149)); */

        /*@ assert (\is_finite((double) (x9))); */
        /*@ assert (x9 < -1.09476442525e-47 || x9 > 1.09476442525e-47); */
        x170 = 1.0/(x9*x9);
        /*@ assert \is_finite((double) (x170)); */
        /*@ assert (x170) >= 0; */
        x171 = 0.25*mu*x144*x170*x71;
        /*@ assert \is_finite((double) (x171)); */
        x172 = 4.0*x10*x9;
        /*@ assert \is_finite((double) (x172)); */
        x173 = -x172*x7;
        /*@ assert \is_finite((double) (x173)); */

        /*@ assert (\is_finite((double) (rt2))); */
        x174 = rt2*rt2*rt2*rt2;
        /*@ assert \is_finite((double) (x174)); */
        x175 = mu*rn;
        /*@ assert \is_finite((double) (x175)); */
        x176 = x175 + x18;
        /*@ assert \is_finite((double) (x176)); */
        x177 = Max(0, x10 + x176);
        /*@ assert \is_finite((double) (x177)); */
        /*@ assert (x177) >= 0; */
        x178 = 2.0*x177;
        /*@ assert \is_finite((double) (x178)); */
        x179 = -x174*x178;
        /*@ assert \is_finite((double) (x179)); */

        /*@ assert (\is_finite((double) (rt1))); */
        x180 = rt1*rt1*rt1*rt1*rt1;
        /*@ assert \is_finite((double) (x180)); */
        x181 = 1.414213562373095145474621858738828450441360473632813*x180;
        /*@ assert \is_finite((double) (x181)); */
        x182 = x148*x181;
        /*@ assert \is_finite((double) (x182)); */
        x183 = -x149*x181;
        /*@ assert \is_finite((double) (x183)); */
        x184 = -x10;
        /*@ assert \is_finite((double) (x184)); */
        x185 = Max(0, x176 + x184);
        /*@ assert \is_finite((double) (x185)); */
        /*@ assert (x185) >= 0; */
        x186 = x174*x185;
        /*@ assert \is_finite((double) (x186)); */
        x187 = 2.0*x186;
        /*@ assert \is_finite((double) (x187)); */
        x188 = 1.414213562373095145474621858738828450441360473632813*rt1*x174;
        /*@ assert \is_finite((double) (x188)); */
        x189 = x148*x188;
        /*@ assert \is_finite((double) (x189)); */
        x190 = -x149*x188;
        /*@ assert \is_finite((double) (x190)); */
        x191 = x177*x7*x8;
        /*@ assert \is_finite((double) (x191)); */
        x192 = -2.0*x191;
        /*@ assert \is_finite((double) (x192)); */
        x193 = x7*x8;
        /*@ assert \is_finite((double) (x193)); */
        x194 = 2.0*x185;
        /*@ assert \is_finite((double) (x194)); */
        x195 = x193*x194;
        /*@ assert \is_finite((double) (x195)); */

        /*@ assert (\is_finite((double) (rt1))); */
        x196 = rt1*rt1*rt1;
        /*@ assert \is_finite((double) (x196)); */
        x197 = 2.828427124746190290949243717477656900882720947265625*x196*x8;
        /*@ assert \is_finite((double) (x197)); */
        x198 = x148*x197;
        /*@ assert \is_finite((double) (x198)); */
        x199 = -x149*x197;
        /*@ assert \is_finite((double) (x199)); */
        x200 = x10*x9;
        /*@ assert \is_finite((double) (x200)); */
        x201 = 2.0*x200;
        /*@ assert \is_finite((double) (x201)); */
        x202 = x201*x7;
        /*@ assert \is_finite((double) (x202)); */
        x203 = x148*x202;
        /*@ assert \is_finite((double) (x203)); */
        x204 = x149*x202;
        /*@ assert \is_finite((double) (x204)); */
        x205 = x173 + x179 + x182 + x183 + x187 + x189 + x190 + x192 + x195 + x198 + x199 + x203 + x204;
        /*@ assert \is_finite((double) (x205)); */
        x210 = x185*x7;
        /*@ assert \is_finite((double) (x210)); */
        x235 = -rt2*x172;
        /*@ assert \is_finite((double) (x235)); */

        /*@ assert (\is_finite((double) (rt2))); */
        x236 = rt2*rt2*rt2;
        /*@ assert \is_finite((double) (x236)); */
        x237 = x178*x236;
        /*@ assert \is_finite((double) (x237)); */
        x238 = x177*x7;
        /*@ assert \is_finite((double) (x238)); */
        x239 = x238*x98;
        /*@ assert \is_finite((double) (x239)); */

        /*@ assert (\is_finite((double) (rt1))); */
        x240 = rt1*rt1*rt1*rt1;
        /*@ assert \is_finite((double) (x240)); */
        x241 = 1.414213562373095145474621858738828450441360473632813*x240;
        /*@ assert \is_finite((double) (x241)); */
        x242 = x148*x241;
        /*@ assert \is_finite((double) (x242)); */
        x243 = -x149*x241;
        /*@ assert \is_finite((double) (x243)); */
        x244 = -x194*x236;
        /*@ assert \is_finite((double) (x244)); */
        x245 = 1.414213562373095145474621858738828450441360473632813*x174;
        /*@ assert \is_finite((double) (x245)); */
        x246 = x148*x245;
        /*@ assert \is_finite((double) (x246)); */
        x247 = -x149*x245;
        /*@ assert \is_finite((double) (x247)); */
        x248 = -x210*x98;
        /*@ assert \is_finite((double) (x248)); */
        x249 = 2.828427124746190290949243717477656900882720947265625*x7*x8;
        /*@ assert \is_finite((double) (x249)); */
        x250 = x148*x249;
        /*@ assert \is_finite((double) (x250)); */
        x251 = -x149*x249;
        /*@ assert \is_finite((double) (x251)); */
        x252 = x148*x201;
        /*@ assert \is_finite((double) (x252)); */
        x253 = rt2*x252;
        /*@ assert \is_finite((double) (x253)); */
        x254 = rt2*x149;
        /*@ assert \is_finite((double) (x254)); */
        x255 = x201*x254;
        /*@ assert \is_finite((double) (x255)); */
        x256 = x235 + x237 + x239 + x242 + x243 + x244 + x246 + x247 + x248 + x250 + x251 + x253 + x255;
        /*@ assert \is_finite((double) (x256)); */
        x295 = x200*x7;
        /*@ assert \is_finite((double) (x295)); */
        x296 = x200*x8;
        /*@ assert \is_finite((double) (x296)); */
        x304 = x10*x148*x7*x9;
        /*@ assert \is_finite((double) (x304)); */
        x356 = 0.25*mu*x170*x71;
        /*@ assert \is_finite((double) (x356)); */
        x357 = 4.0*x144;
        /*@ assert \is_finite((double) (x357)); */
        x358 = x180*x357;
        /*@ assert \is_finite((double) (x358)); */
        x359 = rt1*x174*x357;
        /*@ assert \is_finite((double) (x359)); */
        x360 = x196*x8;
        /*@ assert \is_finite((double) (x360)); */
        x361 = 8.0*x144;
        /*@ assert \is_finite((double) (x361)); */
        x362 = x360*x361;
        /*@ assert \is_finite((double) (x362)); */
        x363 = 2.0*x180;
        /*@ assert \is_finite((double) (x363)); */
        x364 = -x32*x363;
        /*@ assert \is_finite((double) (x364)); */
        x365 = -x35*x363;
        /*@ assert \is_finite((double) (x365)); */
        x366 = 2.0*rt1*x174;
        /*@ assert \is_finite((double) (x366)); */
        x367 = -x32*x366;
        /*@ assert \is_finite((double) (x367)); */
        x368 = -x35*x366;
        /*@ assert \is_finite((double) (x368)); */
        x369 = -4.0*x32*x360;
        /*@ assert \is_finite((double) (x369)); */
        x370 = 4.0*x35;
        /*@ assert \is_finite((double) (x370)); */
        x371 = -x360*x370;
        /*@ assert \is_finite((double) (x371)); */
        x372 = rt2*x10*x9;
        /*@ assert \is_finite((double) (x372)); */
        x373 = -4.0*rt1*x144*x372;
        /*@ assert \is_finite((double) (x373)); */
        x374 = 2.0*x144*x180;
        /*@ assert \is_finite((double) (x374)); */
        x375 = -x148*x374;
        /*@ assert \is_finite((double) (x375)); */
        x376 = -x149*x374;
        /*@ assert \is_finite((double) (x376)); */

        /*@ assert (\is_finite((double) (rt2))); */
        x377 = rt2*rt2*rt2*rt2*rt2;
        /*@ assert \is_finite((double) (x377)); */
        x378 = x144*x148*x377;
        /*@ assert \is_finite((double) (x378)); */
        x379 = 1.414213562373095145474621858738828450441360473632813*x378;
        /*@ assert \is_finite((double) (x379)); */
        x380 = x144*x149*x377;
        /*@ assert \is_finite((double) (x380)); */
        x381 = -1.414213562373095145474621858738828450441360473632813*x380;
        /*@ assert \is_finite((double) (x381)); */
        x382 = 2.0*x144;
        /*@ assert \is_finite((double) (x382)); */
        x383 = x174*x382;
        /*@ assert \is_finite((double) (x383)); */
        x384 = x148*x383;
        /*@ assert \is_finite((double) (x384)); */
        x385 = -rt1*x384;
        /*@ assert \is_finite((double) (x385)); */
        x386 = rt1*x149;
        /*@ assert \is_finite((double) (x386)); */
        x387 = -x383*x386;
        /*@ assert \is_finite((double) (x387)); */
        x388 = rt2*x144*x148*x240;
        /*@ assert \is_finite((double) (x388)); */
        x389 = 1.414213562373095145474621858738828450441360473632813*x388;
        /*@ assert \is_finite((double) (x389)); */
        x390 = rt2*x144*x149*x240;
        /*@ assert \is_finite((double) (x390)); */
        x391 = -1.414213562373095145474621858738828450441360473632813*x390;
        /*@ assert \is_finite((double) (x391)); */
        x392 = 2.0*rt1*x144*x236;
        /*@ assert \is_finite((double) (x392)); */
        x393 = Max(0, x31);
        /*@ assert \is_finite((double) (x393)); */
        /*@ assert (x393) >= 0; */
        x394 = -x392*x393;
        /*@ assert \is_finite((double) (x394)); */
        x395 = Max(0, x34);
        /*@ assert \is_finite((double) (x395)); */
        /*@ assert (x395) >= 0; */
        x396 = x392*x395;
        /*@ assert \is_finite((double) (x396)); */
        x397 = 2.0*rt2*x144*x196;
        /*@ assert \is_finite((double) (x397)); */
        x398 = -x393*x397;
        /*@ assert \is_finite((double) (x398)); */
        x399 = x395*x397;
        /*@ assert \is_finite((double) (x399)); */
        x400 = x144*x148*x236*x7;
        /*@ assert \is_finite((double) (x400)); */
        x401 = 2.828427124746190290949243717477656900882720947265625*x400;
        /*@ assert \is_finite((double) (x401)); */
        x402 = x144*x149*x236*x7;
        /*@ assert \is_finite((double) (x402)); */
        x403 = -2.828427124746190290949243717477656900882720947265625*x402;
        /*@ assert \is_finite((double) (x403)); */
        x404 = 4.0*x144*x196*x8;
        /*@ assert \is_finite((double) (x404)); */
        x405 = -x148*x404;
        /*@ assert \is_finite((double) (x405)); */
        x406 = -x149*x404;
        /*@ assert \is_finite((double) (x406)); */
        x407 = x295*x76;
        /*@ assert \is_finite((double) (x407)); */
        x408 = 1.414213562373095145474621858738828450441360473632813*x35;
        /*@ assert \is_finite((double) (x408)); */
        x409 = x295*x408;
        /*@ assert \is_finite((double) (x409)); */
        x410 = -x409;
        /*@ assert \is_finite((double) (x410)); */
        x411 = x296*x76;
        /*@ assert \is_finite((double) (x411)); */
        x412 = x296*x408;
        /*@ assert \is_finite((double) (x412)); */
        x413 = -x412;
        /*@ assert \is_finite((double) (x413)); */
        x414 = x372*x382;
        /*@ assert \is_finite((double) (x414)); */
        x415 = x148*x414;
        /*@ assert \is_finite((double) (x415)); */
        x416 = rt1*x415;
        /*@ assert \is_finite((double) (x416)); */
        x417 = x386*x414;
        /*@ assert \is_finite((double) (x417)); */
        x418 = 1.414213562373095145474621858738828450441360473632813*x144;
        /*@ assert \is_finite((double) (x418)); */
        x419 = x304*x418;
        /*@ assert \is_finite((double) (x419)); */
        x420 = -x419;
        /*@ assert \is_finite((double) (x420)); */
        x421 = 1.414213562373095145474621858738828450441360473632813*x144*x149;
        /*@ assert \is_finite((double) (x421)); */
        x422 = x295*x421;
        /*@ assert \is_finite((double) (x422)); */
        x423 = x10*x148*x8*x9;
        /*@ assert \is_finite((double) (x423)); */
        x424 = -x418*x423;
        /*@ assert \is_finite((double) (x424)); */
        x425 = x296*x421;
        /*@ assert \is_finite((double) (x425)); */
        x426 = x358 + x359 + x362 + x364 + x365 + x367 + x368 + x369 + x371 + x373 + x375 + x376 + x379 + x381 + x385 + x387 + x389 + x391 + x394 + x396 + x398 + x399 + x401 + x403 + x405 + x406 + x407 + x410 + x411 + x413 + x416 + x417 + x420 + x422 + x424 + x425;
        /*@ assert \is_finite((double) (x426)); */
        x428 = x240*x382;
        /*@ assert \is_finite((double) (x428)); */
        x458 = x144*x377;
        /*@ assert \is_finite((double) (x458)); */
        x459 = -4.0*x458;
        /*@ assert \is_finite((double) (x459)); */
        x460 = rt2*x240;
        /*@ assert \is_finite((double) (x460)); */
        x461 = -x357*x460;
        /*@ assert \is_finite((double) (x461)); */
        x462 = x236*x7;
        /*@ assert \is_finite((double) (x462)); */
        x463 = -x361*x462;
        /*@ assert \is_finite((double) (x463)); */
        x464 = 2.0*x377;
        /*@ assert \is_finite((double) (x464)); */
        x465 = x32*x464;
        /*@ assert \is_finite((double) (x465)); */
        x466 = x35*x464;
        /*@ assert \is_finite((double) (x466)); */
        x467 = 2.0*rt2*x240;
        /*@ assert \is_finite((double) (x467)); */
        x468 = x32*x467;
        /*@ assert \is_finite((double) (x468)); */
        x469 = x35*x467;
        /*@ assert \is_finite((double) (x469)); */
        x470 = x236*x32*x7;
        /*@ assert \is_finite((double) (x470)); */
        x471 = 4.0*x470;
        /*@ assert \is_finite((double) (x471)); */
        x472 = x370*x462;
        /*@ assert \is_finite((double) (x472)); */
        x473 = x296*x357;
        /*@ assert \is_finite((double) (x473)); */
        x474 = 0.5857864376269049655476806037768255919218063354492188*x378;
        /*@ assert \is_finite((double) (x474)); */
        x475 = 3.41421356237309492343001693370752036571502685546875*x380;
        /*@ assert \is_finite((double) (x475)); */
        x476 = -x393*x428;
        /*@ assert \is_finite((double) (x476)); */
        x477 = x395*x428;
        /*@ assert \is_finite((double) (x477)); */
        x478 = 0.5857864376269049655476806037768255919218063354492188*x388;
        /*@ assert \is_finite((double) (x478)); */
        x479 = 3.41421356237309492343001693370752036571502685546875*x390;
        /*@ assert \is_finite((double) (x479)); */
        x480 = 1.171572875253809931095361207553651183843612670898438*x400;
        /*@ assert \is_finite((double) (x480)); */
        x481 = 6.8284271247461898468600338674150407314300537109375*x402;
        /*@ assert \is_finite((double) (x481)); */
        x482 = 2.0*x144*x7*x8;
        /*@ assert \is_finite((double) (x482)); */
        x483 = -x393*x482;
        /*@ assert \is_finite((double) (x483)); */
        x484 = x395*x482;
        /*@ assert \is_finite((double) (x484)); */
        x485 = -x407;
        /*@ assert \is_finite((double) (x485)); */
        x486 = -x411;
        /*@ assert \is_finite((double) (x486)); */
        x487 = -x422;
        /*@ assert \is_finite((double) (x487)); */
        x488 = -0.5857864376269049655476806037768255919218063354492188*x144*x423;
        /*@ assert \is_finite((double) (x488)); */
        x489 = x144*x149;
        /*@ assert \is_finite((double) (x489)); */
        x490 = -3.41421356237309492343001693370752036571502685546875*x10*x489*x8*x9;
        /*@ assert \is_finite((double) (x490)); */
        x491 = x409 + x412 + x419 + x459 + x461 + x463 + x465 + x466 + x468 + x469 + x471 + x472 + x473 + x474 + x475 + x476 + x477 + x478 + x479 + x480 + x481 + x483 + x484 + x485 + x486 + x487 + x488 + x490;
        /*@ assert \is_finite((double) (x491)); */

    }
    if (x84)
    {
        x1 = 1.0*un;
        /*@ assert \is_finite((double) (x1)); */
        x2 = -x1;
        /*@ assert \is_finite((double) (x2)); */
        x3 = mu*rn;
        /*@ assert \is_finite((double) (x3)); */
        x4 = 1.0*x3;
        /*@ assert \is_finite((double) (x4)); */
        x5 = x2 + x4;
        /*@ assert \is_finite((double) (x5)); */
        x18 = -un;
        /*@ assert \is_finite((double) (x18)); */
        x29 = 1.0*x10;
        /*@ assert \is_finite((double) (x29)); */
        x30 = -x29;
        /*@ assert \is_finite((double) (x30)); */
        x34 = x29 + x5;
        /*@ assert \is_finite((double) (x34)); */
        x35 = Heaviside(x34);
        /*@ assert \is_finite((double) (x35)); */

        /*@ assert (\is_finite((double) (x10))); */
        /*@ assert (x10 < -1.09476442525e-47 || x10 > 1.09476442525e-47); */
        x71 = 1.0/x10;
        /*@ assert \is_finite((double) (x71)); */
        x83 = 1.0*mu*x35*x71;
        /*@ assert \is_finite((double) (x83)); */
        x144 = Heaviside(x10);
        /*@ assert \is_finite((double) (x144)); */
        x146 = -x4;
        /*@ assert \is_finite((double) (x146)); */
        x147 = x1 + x146;
        /*@ assert \is_finite((double) (x147)); */
        x148 = Heaviside(x147 + x30);
        /*@ assert \is_finite((double) (x148)); */

        /*@ assert (\is_finite((double) (x9))); */
        /*@ assert (x9 < -1.09476442525e-47 || x9 > 1.09476442525e-47); */
        x170 = 1.0/(x9*x9);
        /*@ assert \is_finite((double) (x170)); */
        /*@ assert (x170) >= 0; */

        /*@ assert (\is_finite((double) (rt2))); */
        x174 = rt2*rt2*rt2*rt2;
        /*@ assert \is_finite((double) (x174)); */
        x175 = mu*rn;
        /*@ assert \is_finite((double) (x175)); */
        x176 = x175 + x18;
        /*@ assert \is_finite((double) (x176)); */
        x177 = Max(0, x10 + x176);
        /*@ assert \is_finite((double) (x177)); */
        /*@ assert (x177) >= 0; */
        x184 = -x10;
        /*@ assert \is_finite((double) (x184)); */
        x185 = Max(0, x176 + x184);
        /*@ assert \is_finite((double) (x185)); */
        /*@ assert (x185) >= 0; */
        x186 = x174*x185;
        /*@ assert \is_finite((double) (x186)); */
        x191 = x177*x7*x8;
        /*@ assert \is_finite((double) (x191)); */
        x193 = x7*x8;
        /*@ assert \is_finite((double) (x193)); */
        x200 = x10*x9;
        /*@ assert \is_finite((double) (x200)); */
        x201 = 2.0*x200;
        /*@ assert \is_finite((double) (x201)); */
        x202 = x201*x7;
        /*@ assert \is_finite((double) (x202)); */
        x203 = x148*x202;
        /*@ assert \is_finite((double) (x203)); */
        x206 = -x202;
        /*@ assert \is_finite((double) (x206)); */
        x207 = x174*x177;
        /*@ assert \is_finite((double) (x207)); */
        x208 = -1.0*x207;
        /*@ assert \is_finite((double) (x208)); */
        x209 = -1.0*x191;
        /*@ assert \is_finite((double) (x209)); */
        x210 = x185*x7;
        /*@ assert \is_finite((double) (x210)); */
        x211 = x210*x8;
        /*@ assert \is_finite((double) (x211)); */
        x212 = x186 + x203 + x206 + x208 + x209 + x211;
        /*@ assert \is_finite((double) (x212)); */

        /*@ assert (\is_finite((double) (rt2))); */
        x236 = rt2*rt2*rt2;
        /*@ assert \is_finite((double) (x236)); */
        x238 = x177*x7;
        /*@ assert \is_finite((double) (x238)); */

        /*@ assert (\is_finite((double) (rt1))); */
        x240 = rt1*rt1*rt1*rt1;
        /*@ assert \is_finite((double) (x240)); */
        x252 = x148*x201;
        /*@ assert \is_finite((double) (x252)); */
        x257 = 0.5*rt1*rt2*x144;
        /*@ assert \is_finite((double) (x257)); */
        x258 = -x201;
        /*@ assert \is_finite((double) (x258)); */
        x259 = x177*x8;
        /*@ assert \is_finite((double) (x259)); */
        x260 = 1.0*x185;
        /*@ assert \is_finite((double) (x260)); */
        x261 = x260*x7;
        /*@ assert \is_finite((double) (x261)); */
        x262 = -x261;
        /*@ assert \is_finite((double) (x262)); */
        x263 = -x260*x8;
        /*@ assert \is_finite((double) (x263)); */
        x264 = x238 + x258 + x259 + x262 + x263;
        /*@ assert \is_finite((double) (x264)); */
        x265 = x252 + x264;
        /*@ assert \is_finite((double) (x265)); */
        x296 = x200*x8;
        /*@ assert \is_finite((double) (x296)); */
        x344 = rt2*x144;
        /*@ assert \is_finite((double) (x344)); */
        x370 = 4.0*x35;
        /*@ assert \is_finite((double) (x370)); */
        x372 = rt2*x10*x9;
        /*@ assert \is_finite((double) (x372)); */

        /*@ assert (\is_finite((double) (rt2))); */
        x377 = rt2*rt2*rt2*rt2*rt2;
        /*@ assert \is_finite((double) (x377)); */
        x382 = 2.0*x144;
        /*@ assert \is_finite((double) (x382)); */
        x383 = x174*x382;
        /*@ assert \is_finite((double) (x383)); */
        x384 = x148*x383;
        /*@ assert \is_finite((double) (x384)); */
        x395 = Max(0, x34);
        /*@ assert \is_finite((double) (x395)); */
        /*@ assert (x395) >= 0; */
        x414 = x372*x382;
        /*@ assert \is_finite((double) (x414)); */
        x415 = x148*x414;
        /*@ assert \is_finite((double) (x415)); */
        x427 = 0.5*mu*x170*x71;
        /*@ assert \is_finite((double) (x427)); */
        x428 = x240*x382;
        /*@ assert \is_finite((double) (x428)); */
        x429 = 4.0*x144*x7;
        /*@ assert \is_finite((double) (x429)); */
        x430 = x429*x8;
        /*@ assert \is_finite((double) (x430)); */
        x431 = x240*x35;
        /*@ assert \is_finite((double) (x431)); */
        x432 = -2.0*x431;
        /*@ assert \is_finite((double) (x432)); */
        x433 = x174*x35;
        /*@ assert \is_finite((double) (x433)); */
        x434 = -2.0*x433;
        /*@ assert \is_finite((double) (x434)); */
        x435 = -x193*x370;
        /*@ assert \is_finite((double) (x435)); */
        x436 = -x414;
        /*@ assert \is_finite((double) (x436)); */
        x437 = -x148*x428;
        /*@ assert \is_finite((double) (x437)); */
        x438 = x144*x236;
        /*@ assert \is_finite((double) (x438)); */
        x439 = -x260*x438;
        /*@ assert \is_finite((double) (x439)); */
        x440 = x395*x438;
        /*@ assert \is_finite((double) (x440)); */
        x441 = -x384;
        /*@ assert \is_finite((double) (x441)); */
        x442 = -x261*x344;
        /*@ assert \is_finite((double) (x442)); */
        x443 = rt2*x144*x7;
        /*@ assert \is_finite((double) (x443)); */
        x444 = x395*x443;
        /*@ assert \is_finite((double) (x444)); */
        x445 = -x148*x430;
        /*@ assert \is_finite((double) (x445)); */
        x446 = x383 + x415 + x428 + x430 + x432 + x434 + x435 + x436 + x437 + x439 + x440 + x441 + x442 + x444 + x445;
        /*@ assert \is_finite((double) (x446)); */
        x458 = x144*x377;
        /*@ assert \is_finite((double) (x458)); */
        x462 = x236*x7;
        /*@ assert \is_finite((double) (x462)); */
        x464 = 2.0*x377;
        /*@ assert \is_finite((double) (x464)); */
        x466 = x35*x464;
        /*@ assert \is_finite((double) (x466)); */
        x467 = 2.0*rt2*x240;
        /*@ assert \is_finite((double) (x467)); */
        x469 = x35*x467;
        /*@ assert \is_finite((double) (x469)); */
        x472 = x370*x462;
        /*@ assert \is_finite((double) (x472)); */
        x492 = 2.0*x458;
        /*@ assert \is_finite((double) (x492)); */
        x493 = rt2*x428;
        /*@ assert \is_finite((double) (x493)); */
        x494 = x236*x429;
        /*@ assert \is_finite((double) (x494)); */
        x495 = -x466;
        /*@ assert \is_finite((double) (x495)); */
        x496 = -x469;
        /*@ assert \is_finite((double) (x496)); */
        x497 = -x472;
        /*@ assert \is_finite((double) (x497)); */
        x498 = x296*x382;
        /*@ assert \is_finite((double) (x498)); */
        x499 = -x498;
        /*@ assert \is_finite((double) (x499)); */
        x500 = x144*x240;
        /*@ assert \is_finite((double) (x500)); */
        x501 = x185*x500;
        /*@ assert \is_finite((double) (x501)); */
        x502 = -x148*x492;
        /*@ assert \is_finite((double) (x502)); */
        x503 = 1.0*x395;
        /*@ assert \is_finite((double) (x503)); */
        x504 = -x500*x503;
        /*@ assert \is_finite((double) (x504)); */
        x505 = -x148*x493;
        /*@ assert \is_finite((double) (x505)); */
        x506 = x144*x211;
        /*@ assert \is_finite((double) (x506)); */
        x507 = -x148*x494;
        /*@ assert \is_finite((double) (x507)); */
        x508 = -x144*x503*x7*x8;
        /*@ assert \is_finite((double) (x508)); */
        x509 = x148*x498;
        /*@ assert \is_finite((double) (x509)); */
        x510 = x492 + x493 + x494 + x495 + x496 + x497 + x499 + x501 + x502 + x504 + x505 + x506 + x507 + x508 + x509;
        /*@ assert \is_finite((double) (x510)); */

    }
    if (x142)
    {
        x3 = mu*rn;
        /*@ assert \is_finite((double) (x3)); */
        x18 = -un;
        /*@ assert \is_finite((double) (x18)); */
        x19 = x18 + x3;
        /*@ assert \is_finite((double) (x19)); */
        x20 = Heaviside(x19);
        /*@ assert \is_finite((double) (x20)); */
        x23 = -2*x22;
        /*@ assert \is_finite((double) (x23)); */
        x24 = x19 + x23;
        /*@ assert \is_finite((double) (x24)); */
        x25 = Heaviside(x24);
        /*@ assert \is_finite((double) (x25)); */
        x26 = 0.5*x25;
        /*@ assert \is_finite((double) (x26)); */
        x66 = mu*ut1;
        /*@ assert \is_finite((double) (x66)); */

        /*@ assert (\is_finite((double) (x15))); */
        /*@ assert (x15 < -1.09476442525e-47 || x15 > 1.09476442525e-47); */
        x67 = 1.0/x15;
        /*@ assert \is_finite((double) (x67)); */
        x68 = 1.0*x67;
        /*@ assert \is_finite((double) (x68)); */
        x95 = mu*ut2;
        /*@ assert \is_finite((double) (x95)); */
        x117 = 0.5*x20*x67;
        /*@ assert \is_finite((double) (x117)); */
        x118 = ut1*x117;
        /*@ assert \is_finite((double) (x118)); */
        x119 = 0.5*x25*x67;
        /*@ assert \is_finite((double) (x119)); */
        x120 = ut1*x119;
        /*@ assert \is_finite((double) (x120)); */
        x133 = ut2*x117;
        /*@ assert \is_finite((double) (x133)); */
        x134 = ut2*x119;
        /*@ assert \is_finite((double) (x134)); */
        x141 = -x118;
        /*@ assert \is_finite((double) (x141)); */

        /*@ assert (\is_finite((double) (x14))); */
        /*@ assert (x14 < -1.09476442525e-47 || x14 > 1.09476442525e-47); */
        x155 = 1.0/x14;
        /*@ assert \is_finite((double) (x155)); */
        x156 = 1.0*mu*x155*x25;
        /*@ assert \is_finite((double) (x156)); */
        x157 = x12*x156;
        /*@ assert \is_finite((double) (x157)); */
        x158 = Max(0, x19);
        /*@ assert \is_finite((double) (x158)); */
        /*@ assert (x158) >= 0; */
        x159 = 0.5*x158;
        /*@ assert \is_finite((double) (x159)); */
        x160 = -x68;
        /*@ assert \is_finite((double) (x160)); */
        x161 = 1.0*x155*x67;
        /*@ assert \is_finite((double) (x161)); */
        x162 = x12*x161;
        /*@ assert \is_finite((double) (x162)); */
        x163 = x160 + x162;
        /*@ assert \is_finite((double) (x163)); */
        x164 = -x159*x163;
        /*@ assert \is_finite((double) (x164)); */
        x165 = Max(0, x24);
        /*@ assert \is_finite((double) (x165)); */
        /*@ assert (x165) >= 0; */
        x166 = 0.5*x165;
        /*@ assert \is_finite((double) (x166)); */
        x167 = -x162;
        /*@ assert \is_finite((double) (x167)); */
        x168 = x167 + x68;
        /*@ assert \is_finite((double) (x168)); */
        x169 = -x166*x168;
        /*@ assert \is_finite((double) (x169)); */
        x229 = ut1*ut2*x155;
        /*@ assert \is_finite((double) (x229)); */
        x230 = 0.5*x158*x229*x67;
        /*@ assert \is_finite((double) (x230)); */
        x231 = -x230;
        /*@ assert \is_finite((double) (x231)); */
        x232 = 1.0*mu*ut1*ut2*x155*x25;
        /*@ assert \is_finite((double) (x232)); */
        x233 = 0.5*x165*x229*x67;
        /*@ assert \is_finite((double) (x233)); */
        x234 = x231 + x232 + x233;
        /*@ assert \is_finite((double) (x234)); */
        x274 = x117*x66;
        /*@ assert \is_finite((double) (x274)); */
        x275 = -x119*x66;
        /*@ assert \is_finite((double) (x275)); */
        x280 = 0.5*x155*x20;
        /*@ assert \is_finite((double) (x280)); */
        x281 = -x12*x280;
        /*@ assert \is_finite((double) (x281)); */
        x282 = 0.5*x155*x25;
        /*@ assert \is_finite((double) (x282)); */
        x283 = -x12*x282;
        /*@ assert \is_finite((double) (x283)); */

        /*@ assert (\is_finite((double) (mu))); */
        /*@ assert (mu < -1.09476442525e-47 || mu > 1.09476442525e-47); */
        x284 = 1.0/mu;
        /*@ assert \is_finite((double) (x284)); */
        x285 = x284*x68;
        /*@ assert \is_finite((double) (x285)); */
        x286 = 1.0*x155*x284*x67;
        /*@ assert \is_finite((double) (x286)); */
        x287 = x12*x286;
        /*@ assert \is_finite((double) (x287)); */
        x288 = -x287;
        /*@ assert \is_finite((double) (x288)); */
        x289 = x285 + x288;
        /*@ assert \is_finite((double) (x289)); */
        x290 = -x159*x289;
        /*@ assert \is_finite((double) (x290)); */
        x291 = -x285;
        /*@ assert \is_finite((double) (x291)); */
        x292 = x287 + x291;
        /*@ assert \is_finite((double) (x292)); */
        x293 = -x166*x292;
        /*@ assert \is_finite((double) (x293)); */
        x324 = -ut1*ut2*x280;
        /*@ assert \is_finite((double) (x324)); */
        x325 = -x229*x26;
        /*@ assert \is_finite((double) (x325)); */
        x326 = x230*x284;
        /*@ assert \is_finite((double) (x326)); */
        x327 = -x233*x284;
        /*@ assert \is_finite((double) (x327)); */
        x328 = x324 + x325 + x326 + x327;
        /*@ assert \is_finite((double) (x328)); */
        x338 = -x133;
        /*@ assert \is_finite((double) (x338)); */
        x451 = x13*x156;
        /*@ assert \is_finite((double) (x451)); */
        x452 = x13*x161;
        /*@ assert \is_finite((double) (x452)); */
        x453 = x160 + x452;
        /*@ assert \is_finite((double) (x453)); */
        x454 = -x159*x453;
        /*@ assert \is_finite((double) (x454)); */
        x455 = -x452;
        /*@ assert \is_finite((double) (x455)); */
        x456 = x455 + x68;
        /*@ assert \is_finite((double) (x456)); */
        x457 = -x166*x456;
        /*@ assert \is_finite((double) (x457)); */
        x519 = x117*x95;
        /*@ assert \is_finite((double) (x519)); */
        x520 = -x119*x95;
        /*@ assert \is_finite((double) (x520)); */
        x554 = -x13*x280;
        /*@ assert \is_finite((double) (x554)); */
        x555 = -x13*x282;
        /*@ assert \is_finite((double) (x555)); */
        x556 = x13*x286;
        /*@ assert \is_finite((double) (x556)); */
        x557 = -x556;
        /*@ assert \is_finite((double) (x557)); */
        x558 = x285 + x557;
        /*@ assert \is_finite((double) (x558)); */
        x559 = -x159*x558;
        /*@ assert \is_finite((double) (x559)); */
        x560 = x291 + x556;
        /*@ assert \is_finite((double) (x560)); */
        x561 = -x166*x560;
        /*@ assert \is_finite((double) (x561)); */

    }
    if (x143)
    {
        x3 = mu*rn;
        /*@ assert \is_finite((double) (x3)); */
        x18 = -un;
        /*@ assert \is_finite((double) (x18)); */
        x19 = x18 + x3;
        /*@ assert \is_finite((double) (x19)); */
        x20 = Heaviside(x19);
        /*@ assert \is_finite((double) (x20)); */
        x21 = 0.5*x20;
        /*@ assert \is_finite((double) (x21)); */
        x23 = -2*x22;
        /*@ assert \is_finite((double) (x23)); */
        x24 = x19 + x23;
        /*@ assert \is_finite((double) (x24)); */
        x25 = Heaviside(x24);
        /*@ assert \is_finite((double) (x25)); */
        x26 = 0.5*x25;
        /*@ assert \is_finite((double) (x26)); */
        x66 = mu*ut1;
        /*@ assert \is_finite((double) (x66)); */

        /*@ assert (\is_finite((double) (x15))); */
        /*@ assert (x15 < -1.09476442525e-47 || x15 > 1.09476442525e-47); */
        x67 = 1.0/x15;
        /*@ assert \is_finite((double) (x67)); */
        x68 = 1.0*x67;
        /*@ assert \is_finite((double) (x68)); */
        x69 = x66*x68;
        /*@ assert \is_finite((double) (x69)); */
        x70 = x25*x69;
        /*@ assert \is_finite((double) (x70)); */
        x95 = mu*ut2;
        /*@ assert \is_finite((double) (x95)); */
        x96 = x68*x95;
        /*@ assert \is_finite((double) (x96)); */
        x97 = x25*x96;
        /*@ assert \is_finite((double) (x97)); */
        x110 = -mu*x21;
        /*@ assert \is_finite((double) (x110)); */
        x111 = mu*x26;
        /*@ assert \is_finite((double) (x111)); */
        x117 = 0.5*x20*x67;
        /*@ assert \is_finite((double) (x117)); */
        x118 = ut1*x117;
        /*@ assert \is_finite((double) (x118)); */
        x119 = 0.5*x25*x67;
        /*@ assert \is_finite((double) (x119)); */
        x120 = ut1*x119;
        /*@ assert \is_finite((double) (x120)); */
        x133 = ut2*x117;
        /*@ assert \is_finite((double) (x133)); */
        x134 = ut2*x119;
        /*@ assert \is_finite((double) (x134)); */
        x339 = -x26;
        /*@ assert \is_finite((double) (x339)); */

    }
    /*@ assert x17 || x28 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[0]; */
        result[0] = x6;
        /*@ assert \is_finite((double) (result[0])); */
    }
    else if (x28)
    {
        /*@ assigns result[0]; */
        result[0] = x21 + x26;
        /*@ assert \is_finite((double) (result[0])); */
    }
    else if (x53)
    {
        /*@ assigns result[0]; */
        result[0] = x33 + x36;
        /*@ assert \is_finite((double) (result[0])); */
    }
    else if (x63)
    {
        /*@ assigns result[0]; */
        result[0] = x58 + x62;
        /*@ assert \is_finite((double) (result[0])); */
    }
    /*@ assert \is_finite((double) (result[0])); */

    /*@ assert x17 || x142 || x143 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[1]; */
        result[1] = 0.0;
        /*@ assert \is_finite((double) (result[1])); */
        /*@ assert (result[1]) >= 0; */
    }
    else if (x142)
    {
        /*@ assigns result[1]; */
        result[1] = x120 + x141;
        /*@ assert \is_finite((double) (result[1])); */
    }
    else if (x143)
    {
        /*@ assigns result[1]; */
        result[1] = 0;
        /*@ assert \is_finite((double) (result[1])); */
        /*@ assert (result[1]) >= 0; */
    }
    else if (x53)
    {
        /*@ assigns result[1]; */
        result[1] = -x153;
        /*@ assert \is_finite((double) (result[1])); */
    }
    else if (x63)
    {
        /*@ assigns result[1]; */
        result[1] = x129 + x154;
        /*@ assert \is_finite((double) (result[1])); */
    }
    /*@ assert \is_finite((double) (result[1])); */

    /*@ assert x17 || x142 || x143 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[2]; */
        result[2] = 0.0;
        /*@ assert \is_finite((double) (result[2])); */
        /*@ assert (result[2]) >= 0; */
    }
    else if (x142)
    {
        /*@ assigns result[2]; */
        result[2] = x134 + x338;
        /*@ assert \is_finite((double) (result[2])); */
    }
    else if (x143)
    {
        /*@ assigns result[2]; */
        result[2] = x21 + x339;
        /*@ assert \is_finite((double) (result[2])); */
    }
    else if (x53)
    {
        /*@ assigns result[2]; */
        result[2] = 1.0*x350*x71;
        /*@ assert \is_finite((double) (result[2])); */
    }
    else if (x63)
    {
        /*@ assigns result[2]; */
        result[2] = x138 + x351;
        /*@ assert \is_finite((double) (result[2])); */
    }
    /*@ assert \is_finite((double) (result[2])); */

    /*@ assert x17 || x28 || x82 || x84 || x63; */
    if (x17)
    {
        /*@ assigns result[3]; */
        result[3] = x65;
        /*@ assert \is_finite((double) (result[3])); */
    }
    else if (x28)
    {
        /*@ assigns result[3]; */
        result[3] = x70;
        /*@ assert \is_finite((double) (result[3])); */
    }
    else if (x82)
    {
        /*@ assigns result[3]; */
        result[3] = x72*x81;
        /*@ assert \is_finite((double) (result[3])); */
    }
    else if (x84)
    {
        /*@ assigns result[3]; */
        result[3] = rt1*x83;
        /*@ assert \is_finite((double) (result[3])); */
    }
    else if (x63)
    {
        /*@ assigns result[3]; */
        result[3] = x91 + x94;
        /*@ assert \is_finite((double) (result[3])); */
    }
    /*@ assert \is_finite((double) (result[3])); */

    /*@ assert x17 || x142 || x143 || x82 || x84 || x63; */
    if (x17)
    {
        /*@ assigns result[4]; */
        result[4] = 0.0;
        /*@ assert \is_finite((double) (result[4])); */
        /*@ assert (result[4]) >= 0; */
    }
    else if (x142)
    {
        /*@ assigns result[4]; */
        result[4] = x157 + x164 + x169;
        /*@ assert \is_finite((double) (result[4])); */
    }
    else if (x143)
    {
        /*@ assigns result[4]; */
        result[4] = 0;
        /*@ assert \is_finite((double) (result[4])); */
        /*@ assert (result[4]) >= 0; */
    }
    else if (x82)
    {
        /*@ assigns result[4]; */
        result[4] = -x171*x205;
        /*@ assert \is_finite((double) (result[4])); */
    }
    else if (x84)
    {
        /*@ assigns result[4]; */
        result[4] = -0.5*mu*x144*x170*x212*x71;
        /*@ assert \is_finite((double) (result[4])); */
    }
    else if (x63)
    {
        /*@ assigns result[4]; */
        result[4] = x213 + x214 + x221 + x228;
        /*@ assert \is_finite((double) (result[4])); */
    }
    /*@ assert \is_finite((double) (result[4])); */

    /*@ assert x17 || x142 || x143 || x82 || x84 || x63; */
    if (x17)
    {
        /*@ assigns result[5]; */
        result[5] = x355;
        /*@ assert \is_finite((double) (result[5])); */
    }
    else if (x142)
    {
        /*@ assigns result[5]; */
        result[5] = x234;
        /*@ assert \is_finite((double) (result[5])); */
    }
    else if (x143)
    {
        /*@ assigns result[5]; */
        result[5] = -x70;
        /*@ assert \is_finite((double) (result[5])); */
    }
    else if (x82)
    {
        /*@ assigns result[5]; */
        result[5] = -x356*x426;
        /*@ assert \is_finite((double) (result[5])); */
    }
    else if (x84)
    {
        /*@ assigns result[5]; */
        result[5] = x40*x427*x446;
        /*@ assert \is_finite((double) (result[5])); */
    }
    else if (x63)
    {
        /*@ assigns result[5]; */
        result[5] = x271 + x448 + x449 + x450;
        /*@ assert \is_finite((double) (result[5])); */
    }
    /*@ assert \is_finite((double) (result[5])); */

    /*@ assert x17 || x28 || x82 || x84 || x63; */
    if (x17)
    {
        /*@ assigns result[6]; */
        result[6] = x65;
        /*@ assert \is_finite((double) (result[6])); */
    }
    else if (x28)
    {
        /*@ assigns result[6]; */
        result[6] = x97;
        /*@ assert \is_finite((double) (result[6])); */
    }
    else if (x82)
    {
        /*@ assigns result[6]; */
        result[6] = x101*x72;
        /*@ assert \is_finite((double) (result[6])); */
    }
    else if (x84)
    {
        /*@ assigns result[6]; */
        result[6] = rt2*x83;
        /*@ assert \is_finite((double) (result[6])); */
    }
    else if (x63)
    {
        /*@ assigns result[6]; */
        result[6] = x105 + x108;
        /*@ assert \is_finite((double) (result[6])); */
    }
    /*@ assert \is_finite((double) (result[6])); */

    /*@ assert x17 || x142 || x143 || x82 || x84 || x63; */
    if (x17)
    {
        /*@ assigns result[7]; */
        result[7] = 0.0;
        /*@ assert \is_finite((double) (result[7])); */
        /*@ assert (result[7]) >= 0; */
    }
    else if (x142)
    {
        /*@ assigns result[7]; */
        result[7] = x234;
        /*@ assert \is_finite((double) (result[7])); */
    }
    else if (x143)
    {
        /*@ assigns result[7]; */
        result[7] = 0;
        /*@ assert \is_finite((double) (result[7])); */
        /*@ assert (result[7]) >= 0; */
    }
    else if (x82)
    {
        /*@ assigns result[7]; */
        result[7] = x171*x256*x40;
        /*@ assert \is_finite((double) (result[7])); */
    }
    else if (x84)
    {
        /*@ assigns result[7]; */
        result[7] = -mu*x170*x257*x265*x71;
        /*@ assert \is_finite((double) (result[7])); */
    }
    else if (x63)
    {
        /*@ assigns result[7]; */
        result[7] = x269 + x271 + x272 + x273;
        /*@ assert \is_finite((double) (result[7])); */
    }
    /*@ assert \is_finite((double) (result[7])); */

    /*@ assert x17 || x142 || x143 || x82 || x84 || x63; */
    if (x17)
    {
        /*@ assigns result[8]; */
        result[8] = x355;
        /*@ assert \is_finite((double) (result[8])); */
    }
    else if (x142)
    {
        /*@ assigns result[8]; */
        result[8] = x451 + x454 + x457;
        /*@ assert \is_finite((double) (result[8])); */
    }
    else if (x143)
    {
        /*@ assigns result[8]; */
        result[8] = -x97;
        /*@ assert \is_finite((double) (result[8])); */
    }
    else if (x82)
    {
        /*@ assigns result[8]; */
        result[8] = x356*x491;
        /*@ assert \is_finite((double) (result[8])); */
    }
    else if (x84)
    {
        /*@ assigns result[8]; */
        result[8] = -x427*x510;
        /*@ assert \is_finite((double) (result[8])); */
    }
    else if (x63)
    {
        /*@ assigns result[8]; */
        result[8] = x511 + x512 + x515 + x518;
        /*@ assert \is_finite((double) (result[8])); */
    }
    /*@ assert \is_finite((double) (result[8])); */

    /*@ assert x17 || x28 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[9]; */
        result[9] = mu + x109;
        /*@ assert \is_finite((double) (result[9])); */
    }
    else if (x28)
    {
        /*@ assigns result[9]; */
        result[9] = mu + x110 + x112;
        /*@ assert \is_finite((double) (result[9])); */
    }
    else if (x53)
    {
        /*@ assigns result[9]; */
        result[9] = mu + x113 + x114;
        /*@ assert \is_finite((double) (result[9])); */
    }
    else if (x63)
    {
        /*@ assigns result[9]; */
        result[9] = mu + x115 + x116;
        /*@ assert \is_finite((double) (result[9])); */
    }
    /*@ assert \is_finite((double) (result[9])); */

    /*@ assert x17 || x142 || x143 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[10]; */
        result[10] = 0.0;
        /*@ assert \is_finite((double) (result[10])); */
        /*@ assert (result[10]) >= 0; */
    }
    else if (x142)
    {
        /*@ assigns result[10]; */
        result[10] = x274 + x275;
        /*@ assert \is_finite((double) (result[10])); */
    }
    else if (x143)
    {
        /*@ assigns result[10]; */
        result[10] = 0;
        /*@ assert \is_finite((double) (result[10])); */
        /*@ assert (result[10]) >= 0; */
    }
    else if (x53)
    {
        /*@ assigns result[10]; */
        result[10] = mu*x153;
        /*@ assert \is_finite((double) (result[10])); */
    }
    else if (x63)
    {
        /*@ assigns result[10]; */
        result[10] = x277 + x279;
        /*@ assert \is_finite((double) (result[10])); */
    }
    /*@ assert \is_finite((double) (result[10])); */

    /*@ assert x17 || x142 || x143 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[11]; */
        result[11] = 0.0;
        /*@ assert \is_finite((double) (result[11])); */
        /*@ assert (result[11]) >= 0; */
    }
    else if (x142)
    {
        /*@ assigns result[11]; */
        result[11] = x519 + x520;
        /*@ assert \is_finite((double) (result[11])); */
    }
    else if (x143)
    {
        /*@ assigns result[11]; */
        result[11] = x110 + x111;
        /*@ assert \is_finite((double) (result[11])); */
    }
    else if (x53)
    {
        /*@ assigns result[11]; */
        result[11] = -mu*x521*x526;
        /*@ assert \is_finite((double) (result[11])); */
    }
    else if (x63)
    {
        /*@ assigns result[11]; */
        result[11] = x527 + x528;
        /*@ assert \is_finite((double) (result[11])); */
    }
    /*@ assert \is_finite((double) (result[11])); */

    /*@ assert x17 || x28 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[12]; */
        result[12] = 0.0;
        /*@ assert \is_finite((double) (result[12])); */
        /*@ assert (result[12]) >= 0; */
    }
    else if (x28)
    {
        /*@ assigns result[12]; */
        result[12] = x118 + x121;
        /*@ assert \is_finite((double) (result[12])); */
    }
    else if (x53)
    {
        /*@ assigns result[12]; */
        result[12] = rt1*x125;
        /*@ assert \is_finite((double) (result[12])); */
    }
    else if (x63)
    {
        /*@ assigns result[12]; */
        result[12] = x130 + x132;
        /*@ assert \is_finite((double) (result[12])); */
    }
    /*@ assert \is_finite((double) (result[12])); */

    /*@ assert x17 || x142 || x143 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[13]; */
        result[13] = 1.0;
        /*@ assert \is_finite((double) (result[13])); */
        /*@ assert (result[13]) >= 0; */
        /*@ assert (result[13]) != 0; */
    }
    else if (x142)
    {
        /*@ assigns result[13]; */
        result[13] = x281 + x283 + x290 + x293 + 1;
        /*@ assert \is_finite((double) (result[13])); */
    }
    else if (x143)
    {
        /*@ assigns result[13]; */
        result[13] = 1;
        /*@ assert \is_finite((double) (result[13])); */
        /*@ assert (result[13]) >= 0; */
        /*@ assert (result[13]) != 0; */
    }
    else if (x53)
    {
        /*@ assigns result[13]; */
        result[13] = x294*x308;
        /*@ assert \is_finite((double) (result[13])); */
    }
    else if (x63)
    {
        /*@ assigns result[13]; */
        result[13] = x311 + x314 + x320 + x323 + 1;
        /*@ assert \is_finite((double) (result[13])); */
    }
    /*@ assert \is_finite((double) (result[13])); */

    /*@ assert x17 || x142 || x143 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[14]; */
        result[14] = x529 + 0.7071067811865475727373109293694142252206802368164063;
        /*@ assert \is_finite((double) (result[14])); */
    }
    else if (x142)
    {
        /*@ assigns result[14]; */
        result[14] = x328;
        /*@ assert \is_finite((double) (result[14])); */
    }
    else if (x143)
    {
        /*@ assigns result[14]; */
        result[14] = x118 + x120;
        /*@ assert \is_finite((double) (result[14])); */
    }
    else if (x53)
    {
        /*@ assigns result[14]; */
        result[14] = x170*x40*x521*x551;
        /*@ assert \is_finite((double) (result[14])); */
    }
    else if (x63)
    {
        /*@ assigns result[14]; */
        result[14] = x334 + x552 + x553;
        /*@ assert \is_finite((double) (result[14])); */
    }
    /*@ assert \is_finite((double) (result[14])); */

    /*@ assert x17 || x28 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[15]; */
        result[15] = 0.0;
        /*@ assert \is_finite((double) (result[15])); */
        /*@ assert (result[15]) >= 0; */
    }
    else if (x28)
    {
        /*@ assigns result[15]; */
        result[15] = x133 + x135;
        /*@ assert \is_finite((double) (result[15])); */
    }
    else if (x53)
    {
        /*@ assigns result[15]; */
        result[15] = rt2*x125;
        /*@ assert \is_finite((double) (result[15])); */
    }
    else if (x63)
    {
        /*@ assigns result[15]; */
        result[15] = x139 + x140;
        /*@ assert \is_finite((double) (result[15])); */
    }
    /*@ assert \is_finite((double) (result[15])); */

    /*@ assert x17 || x142 || x143 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[16]; */
        result[16] = 0.0;
        /*@ assert \is_finite((double) (result[16])); */
        /*@ assert (result[16]) >= 0; */
    }
    else if (x142)
    {
        /*@ assigns result[16]; */
        result[16] = x328;
        /*@ assert \is_finite((double) (result[16])); */
    }
    else if (x143)
    {
        /*@ assigns result[16]; */
        result[16] = 0;
        /*@ assert \is_finite((double) (result[16])); */
        /*@ assert (result[16]) >= 0; */
    }
    else if (x53)
    {
        /*@ assigns result[16]; */
        result[16] = x170*x257*x331*x71;
        /*@ assert \is_finite((double) (result[16])); */
    }
    else if (x63)
    {
        /*@ assigns result[16]; */
        result[16] = x334 + x336 + x337;
        /*@ assert \is_finite((double) (result[16])); */
    }
    /*@ assert \is_finite((double) (result[16])); */

    /*@ assert x17 || x142 || x143 || x53 || x63; */
    if (x17)
    {
        /*@ assigns result[17]; */
        result[17] = x529 + 1.707106781186547461715008466853760182857513427734375;
        /*@ assert \is_finite((double) (result[17])); */
    }
    else if (x142)
    {
        /*@ assigns result[17]; */
        result[17] = x554 + x555 + x559 + x561 + 1;
        /*@ assert \is_finite((double) (result[17])); */
    }
    else if (x143)
    {
        /*@ assigns result[17]; */
        result[17] = x133 + x134 + 1;
        /*@ assert \is_finite((double) (result[17])); */
    }
    else if (x53)
    {
        /*@ assigns result[17]; */
        result[17] = x294*x585;
        /*@ assert \is_finite((double) (result[17])); */
    }
    else if (x63)
    {
        /*@ assigns result[17]; */
        result[17] = x587 + x589 + x592 + x595 + 1;
        /*@ assert \is_finite((double) (result[17])); */
    }
    /*@ assert \is_finite((double) (result[17])); */
}
#ifdef __FRAMAC__
int main()
{
    double rn =  Frama_C_double_interval(0.0, 1.0e+6);
    double rt1 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double rt2 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double un =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double ut1 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double ut2 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double mu =  Frama_C_double_interval(0.0, 1.0);
    double rhon =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double rhot1 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double rhot2 =  Frama_C_double_interval(-1.0e+6, 1.0e+6);
    double result[18];
    fc3d_NaturalMapABGenerated(rn, rt1, rt2, un, ut1, ut2, mu, rhon, rhot1, rhot2, result);
}
#endif
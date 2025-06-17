// -*- C++ -*-

#ifndef TPC_PAD_HELPER_HH
#define TPC_PAD_HELPER_HH

#include <iostream>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <sstream>
#include <iomanip>

#include <TDirectory.h>
#include <TH2Poly.h>
#include <TMath.h>
#include <TVector3.h>

#include "DetectorID.hh"

namespace tpc
{
const Double_t ZTarget = -143.; // Target from center
  //const Double_t TargetVtxWindow = 25.;
const Double_t TargetVtxWindow = 30.;

enum EPadParameter
{
  kLayerID,
  kNumOfPad,
  kRadius,
  kNumOfDivision,
  kDummy,
  kLength,
  NPadParameter
};

//for kinematic fitting
static const Double_t PositionScale = 1.;

//for clustering
static const Int_t MaxRowDifTPC = 2;
//static const Int_t MaxRowDifTPC = 3;
//for helix tracking
static const Double_t ConstC = 0.299792458; //=c/10^9
//_____________________________________________________________________________
//#OfPad #division #radius padLength
static const Double_t padParameter[NumOfLayersTPC][NPadParameter] =
{{0, 48,    14.75, 48, 0,  9.},
 {1, 48,    24.25, 48, 0,  9.},
 {2, 72,    33.75, 72, 0,  9.},
 {3, 96,    43.25, 96, 0,  9.},
 {4, 120,    52.75,120,0,   9.},
 {5, 144,    62.25,144,0,   9.},
 {6, 168,    71.75,168,0,   9.},
 {7, 192,    81.25,192,0,   9.},
 {8, 216,    90.75,216,0,   9.},
 {9, 240,    100.25,240,0,  9.},
 {10,208,    111.5,241, 0,  12.5},
 {11,218,    124.5,271, 0,  12.5},
 {12,230,    137.5,300, 0,  12.5},
 {13,214,    150.5,330, 0,  12.5},
 {14,212,    163.5,360, 0,  12.5},
 {15,214,    176.5,390, 0,  12.5},
 {16,220,    189.5,420, 0,  12.5},
 {17,224,    202.5,449, 0,  12.5},
 {18,232,    215.5,479, 0,  12.5},
 {19,238,    228.5,509, 0,  12.5},
 {20,244,    241.5,539, 0,  12.5},
 {21,232,    254.5,569, 0,  12.5},
 {22,218,    267.5,599, 0,  12.5},
 {23,210,    280.5,628, 0,  12.5},
 {24,206,    293.5,658, 0,  12.5},
 {25,202,    306.5,688, 0,  12.5},
 {26,200,    319.5,718, 0,  12.5},
 {27,196,    332.5,748, 0,  12.5},
 {28,178,    345.5,777, 0,  12.5},
 {29,130,    358.5,807, 0,  12.5},
 {30,108,    371.5,837, 0,  12.5},
 {31,90,     384.5,867, 0, 12.5}};

//_____________________________________________________________________________
static const Int_t noisypad[] =
{
  1507,4106,4356,4592,5113,5216
};

//_____________________________________________________________________________
static const Int_t padOnCenterFrame[] =
{

//Pads on the frame
965,966,967,968,969,970,971,972,973,974,1019,1020,1021,1022,1023,1024,1025,1026,1027,1028,1177,1178,1179,1180,1181,1182,1183,1184,1185,1186,1187,1188,1189,1190,1191,1192,1193,1194,1195,1196,1197,1198,1199,1200,1201,1202,1203,1204,1205,1206,1207,1208,1209,1210,1211,1212,1237,1238,1239,1240,1241,1242,1243,1244,1245,1246,1247,1248,1249,1250,1251,1252,1253,1254,1255,1256,1257,1258,1259,1260,1261,1262,1263,1264,1265,1266,1267,1268,1269,1270,1271,1272,1395,1396,1397,1398,1399,1400,1401,1405,1406,1407,1408,1409,1410,1411,1412,1413,1414,1415,1416,1417,1418,1419,1420,1421,1422,1423,1424,1425,1426,1427,1428,1429,1430,1431,1432,1436,1437,1438,1439,1440,1441,1442,1455,1456,1457,1458,1459,1460,1461,1465,1466,1467,1468,1469,1470,1471,1472,1473,1474,1475,1476,1477,1478,1479,1480,1481,1482,1483,1484,1485,1486,1487,1488,1489,1490,1491,1492,1496,1497,1498,1499,1500,1501,1502,1595,1596,1597,1598,1604,1605,1606,1607,1608,1647,1648,1649,1650,1651,1652,1658,1659,1660,1663,1664,1665,1671,1672,1673,1674,1675,1676,1715,1716,1717,1718,1719,1725,1726,1727,1728,1807,1808,1809,1810,1811,1812,1813,1814,1815,1816,1817,1818,1879,1880,1881,1882,1889,1890,1891,1892,1953,1954,1955,1956,1957,1958,1959,1960,1961,1962,1963,1964,2018,2019,2020,2025,2026,2027,2028,2100,2101,2102,2103,2104,2111,2112,2113,2114,2115,2187,2188,2189,2190,2195,2196,2197,2219,2220,2221,2222,2223,2224,2225,2226,2227,2228,2309,2310,2311,2317,2318,2323,2324,2330,2331,2332,2413,2414,2415,2416,2417,2418,2419,2420,2421,2422,2427,2428,2429,2430,2518,2519,2520,2521,2522,2523,2524,2525,2526,2527,2540,2541,2542,2543,2544,2545,2546,2547,2548,2549,2637,2638,2639,2640,2732,2733,2739,2740,2761,2762,2768,2769,2950,2951,2952,2953,2954,2955,2956,2957,2958,2987,2988,2989,2990,2991,2992,2993,2994,2995,3174,3175,3176,3177,3178,3179,3180,3181,3182,3219,3220,3221,3222,3223,3224,3225,3226,3227,3405,3406,3407,3408,3409,3410,3411,3412,3413,3458,3459,3460,3461,3462,3463,3464,3465,3466,3642,3643,3644,3645,3646,3647,3648,3649,3650,3703,3704,3705,3706,3707,3708,3709,3710,3711,3877,3878,3879,3880,3881,3882,3883,3884,3945,3946,3947,3948,3949,3950,3951,3952,4098,4099,4105,4174,4180,4181,4308,4309,4310,4311,4312,4313,4314,4315,4316,4391,4392,4393,4394,4395,4396,4397,4398,4399,4512,4513,4520,4603,4610,4611,4713,4714,4715,4716,4717,4718,4719,4720,4811,4812,4813,4814,4815,4816,4817,4818,4910,4911,4912,4913,4914,4915,4916,4917,5016,5017,5018,5019,5020,5021,5022,5023,5104,5105,5108,5111,5112,5217,5218,5221,5224,5225,5287,5288,5289,5290,5291,5292,5293,5294,5295,5408,5409,5410,5411,5412,5413,5414,5415,5416,5441,5442,5443,5444,5445,5566,5567,5568,5569,5570,

//Empty Pads
1017,1018,
1394,1402,1403,1404,1433,1434,1435,1462,1463,1464,1493,1494,1495,1503,1599,1600,1601,1602,1603,1604,1653,1654,1655,1656,1657,1666,1667,1668,1669,1670,1720,1721,1722,1723,1724,1877,1878,1883,1884,1885,1886,1887,1888,1893,1894,2021,2022,2023,2024,2105,2106,2107,2108,2109,2110,2191,2192,2193,2194,2312,2313,2314,2315,2316,2325,2326,2327,2328,2329,2734,2735,2736,2737,2738,2763,2764,2765,2766,2767,4100,4101,4102,4103,4104,4175,4176,4177,4178,4179,4514,4515,4516,4517,4518,4519,4604,4605,4606,4607,4608,4609,5106,5107,5109,5110,5219,5220,5222,5223
};

//_____________________________________________________________________________
static const Int_t deadChannel[NumOfLayersTPC][5] =
{
  //18,179,333,3809, //Dead pads
  //2748, //Empty pads
  //Pads on the gem supporing frame
  //73,135,136,223,333,409,468,627,810,60,111,186,285,556,727,922,1141,1364,1565,1566,1779,1819,1820,2037,2038,2246,2247,2248,2456,2457,2669,2670,2888,3112,3113,3343,3344,3580,3581,3815,4036,1630,1851,2069,2278,2487,2488,2700,2701,2919,3143,3375,3612,3846,4067,4278,4482,4682,4880,5074
  {18, -1, -1, -1, -1}, //0
  {60, 73, -1, -1, -1},
  {111, 135, 136, -1, -1},
  {179, 186, 223, -1, -1},
  {285, 333, -1, -1, -1},
  {409, 468, -1, -1, -1}, //5
  {556, 627, -1, -1, -1},
  {727, 810, -1, -1, -1},
  {922, -1, -1, -1, -1},
  {1141, -1, -1, -1, -1},
  {1364, -1, -1, -1, -1}, //10
  {1565, 1566, 1630, -1, -1},
  {1779, 1819, 1820, 1851, -1},
  {2037, 2038, 2069, -1, -1},
  {2246, 2247, 2248, 2278, -1},
  {2456, 2457, 2487, 2488, -1}, //15
  {2669, 2670, 2700, 2701, 2748},
  {2888, 2919, -1, -1, -1},
  {3112, 3113, 3143, -1, -1},
  {3343, 3344, 3375, -1, -1},
  {3580, 3581, 3612, -1, -1}, //20
  {3809, 3815, 3846, -1, -1},
  {4036, 4067, -1, -1, -1},
  {4278, -1, -1, -1, -1},
  {4481, 4482, -1, -1, -1},
  {4682, -1, -1, -1, -1}, //25
  {4880, -1, -1, -1, -1},
  {5074, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1},
  {-1, -1, -1, -1, -1}, //30
  {-1, -1, -1, -1, -1}

};

//_____________________________________________________________________________
static const Int_t FrameLowEdge[NumOfLayersTPC][5] =
{
  //dummy value : -100
  //Layer0~7 : no pad on the frame
  //Layer10~31(Outer layers) : the last pad treated as a pad on the frame

  {-100,-100,-100,-100,-100}, //0
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100}, //5
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {76,128,-100,-100,-100},
  {72,132,-100,-100,-100},
  {49,110,208,-100,-100}, //10
  {42,94,110,162,218},
  {36,106,182,230,-100},
  {17,99,186,214,-100},
  {4,94,108,198,212},
  {91,113,210,-100,-100}, //15
  {91,120,220,-100,-100},
  {89,126,224,-100,-100},
  {89,134,232,-100,-100},
  {88,141,238,-100,-100},
  {87,148,244,-100,-100}, //20
  {78,146,232,-100,-100},
  {67,143,218,-100,-100},
  {59,142,210,-100,-100},
  {53,144,206,-100,-100},
  {48,146,202,-100,-100}, //25
  {43,149,200,-100,-100},
  {37,150,196,-100,-100},
  {24,145,178,-100,-100},
  //{83,-100,-100,-100,-100},
  //{72,-100,-100,-100,-100}, //30
  //{63,-100,-100,-100,-100}
  {125,-100,-100,-100,-100},
  {108,-100,-100,-100,-100}, //30
  {90,-100,-100,-100,-100}
};

//_____________________________________________________________________________
static const Int_t FrameHighEdge[NumOfLayersTPC][5] =
{
  //dummy value : -100
  //Layer0~7 : no pad on the frame
  //Layer10~31(Outer layers) : the first pad treated as a pad on the frame

  {-100,-100,-100,-100,-100}, //0
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100}, //5
  {-100,-100,-100,-100,-100},
  {-100,-100,-100,-100,-100},
  {85,139,-100,-100,-100},
  {107,167,-100,-100,-100},
  {-1,97,158,-100,-100}, //10
  {-1,55,107,123,175},
  {-1,49,123,193,-100},
  {-1,27,114,196,-100},
  {-1,13,103,117,207},
  {3,100,122,-100,-100}, //15
  {-1,99,128,-100,-100},
  {-1,97,134,-100,-100},
  {-1,97,142,-100,-100},
  {-1,96,149,-100,-100},
  {-1,95,156,-100,-100}, //20
  {-1,85,153,-100,-100},
  {-1,74,150,-100,-100},
  {-1,67,150,-100,-100},
  {-1,61,152,-100,-100},
  {-1,55,153,-100,-100}, //25
  {-1,50,156,-100,-100},
  {-1,45,158,-100,-100},
  {-1,32,153,-100,-100},
  {4,-100,-100,-100,-100},
  {-1,-100,-100,-100,-100}, //30
  {-1,-100,-100,-100,-100}
};
static const Double_t ClusterSizeInner[2][10] ={
  {0.12912,0.12912,0.129121,0.164862,0.211372,0.256798,0.33125,0.315186,0.27439,0.303371},//proton
  {0.274143,0.550639,0.683333,0.683333,0.683333,0.683333,0.683333,0.683333,0.683333,0.683333}//pion
  //0.0 - 0.1 , 0.1 - 0.2, ... 0.9 - 1.0 GeV/c
};
static const Double_t ClusterSizeOuter[2][10] ={
  {0.285714,0.285714,0.326683,0.280962,0.243775,0.234132,0.217313,0.256798,0.293478,0.285319},//proton
  {0.257669,0.348881,0.308642,0.308642,0.308642,0.308642,0.308642,0.308642,0.308642,0.308642}//pion
  //0.0 - 0.1 , 0.1 - 0.2, ... 0.9 - 1.0 GeV/c
};
//_____________________________________________________________________________
inline Int_t GetAGETId(Int_t asad, Int_t layer, Int_t row)
{
  Int_t flag=-1;
  switch(asad){
  case 0:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==1 && row==22-ch*2) flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==46-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==28-(ch-23)*2) flag=0;
    if(layer==1 && row==0) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==1 && row==23-ch*2) flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==47-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==29-(ch-23)*2) flag=1;
    if(layer==1 && row==1) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==0 && row==22-ch*2) flag=2;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==34-(ch-13)*2) flag=2;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==16-(ch-23)*2) flag=2;
    if(layer==0 && row==0) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==0 && row==23-ch*2) flag=3;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==35-(ch-13)*2) flag=3;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==17-(ch-23)*2) flag=3;
    if(layer==0 && row==1) flag=3;
    return flag;
  case 1:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==1 && row==47-ch*2) flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==95-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==77-(ch-23)*2) flag=0;
    if(layer==1 && row==25) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==1 && row==46-ch*2) flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==3 && row==94-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==3 && row==76-(ch-23)*2) flag=1;
    if(layer==1 && row==24) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==0 && row==47-ch*2) flag=2;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==71-(ch-13)*2) flag=2;
    for(Int_t ch=23; ch<=37; ++ch) if(layer==2 && row==53-(ch-23)*2) flag=2;
    if(layer==0 && row==25) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==0 && row==46-ch*2) flag=3;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==2 && row==70-(ch-13)*2) flag=3;
    for(Int_t ch=23; ch<=31; ++ch) if(layer==2 && row==52-(ch-23)*2) flag=3;
    if(layer==0 && row==24) flag=3;
    return flag;
  case 2:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==22-ch*2) flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==5 && row==143-(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==125-(ch-23)*2) flag=0;
    if(layer==5 && row==0) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==23-ch*2) flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==5 && row==142-(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==124-(ch-23)*2) flag=1;
    if(layer==5 && row==1) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=9; ++ch) if(layer==4 && row==18-ch*2) flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==117-(ch-12)*2) flag=2;
    if(layer==4 && row==119) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=9; ++ch) if(layer==4 && row==19-ch*2) flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==116-(ch-12)*2) flag=3;
    if(layer==4 && row==118) flag=3;
    return flag;
  case 3:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==24+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==46+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==66+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==25+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==47+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==67+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==4 && row==20+ch*2) flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==42+(ch-12)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==4 && row==21+ch*2) flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==43+(ch-12)*2) flag=3;
    return flag;
  case 4:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==73+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==95+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==115+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==5 && row==72+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==5 && row==94+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==5 && row==114+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==4 && row==61+ch*2) flag=2;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==83+(ch-12)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==4 && row==60+ch*2) flag=3;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==4 && row==82+(ch-12)*2) flag=3;
    return flag;
  case 5:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==30-ch*2) flag=0;
    for(Int_t ch=12; ch<=16; ++ch) if(layer==7 && row==8-(ch-12)*2) flag=0;
    for(Int_t ch=17; ch<=21; ++ch) if(layer==7 && row==191-(ch-17)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==181-(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==31-ch*2) flag=1;
    for(Int_t ch=12; ch<=16; ++ch) if(layer==7 && row==9-(ch-12)*2) flag=1;
    for(Int_t ch=17; ch<=21; ++ch) if(layer==7 && row==190-(ch-17)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==180-(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==26-ch*2) flag=2;
    for(Int_t ch=12; ch<=14; ++ch) if(layer==6 && row==4-(ch-12)*2) flag=2;
    for(Int_t ch=15; ch<=21; ++ch) if(layer==6 && row==167-(ch-15)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==153-(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==27-ch*2) flag=3;
    for(Int_t ch=12; ch<=14; ++ch) if(layer==6 && row==5-(ch-12)*2) flag=3;
    for(Int_t ch=15; ch<=21; ++ch) if(layer==6 && row==166-(ch-15)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==152-(ch-23)*2) flag=3;
    return flag;
  case 6:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==32+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==54+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==74+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==33+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==55+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==75+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==28+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==50+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==70+(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==29+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==51+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==71+(ch-23)*2) flag=3;
    return flag;
  case 7:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==97+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==119+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==139+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==7 && row==96+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==7 && row==118+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==7 && row==138+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==85+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==107+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==127+(ch-23)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==6 && row==84+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==6 && row==106+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==6 && row==126+(ch-23)*2) flag=3;
    return flag;
  case 8:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==9 && row==ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==22+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==42+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==86+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==9 && row==1+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==23+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==43+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==87+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==8 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==86+(ch-46)*2) flag=2;
    if(layer==8 && row==106) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==8 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==87+(ch-46)*2) flag=3;
    if(layer==8 && row==107) flag=3;
    return flag;
  case 9:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==9 && row==121+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==143+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==163+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==207+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==227+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==9 && row==120+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==9 && row==142+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==9 && row==162+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==9 && row==206+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==9 && row==226+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==8 && row==109+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==131+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==151+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==195+(ch-46)*2) flag=2;
    if(layer==8 && row==215) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==8 && row==108+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==8 && row==130+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==8 && row==150+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==8 && row==194+(ch-46)*2) flag=3;
    if(layer==8 && row==214) flag=3;
    return flag;
  case 10:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==11 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==11 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==87+(ch-46)*2) flag=0;
    if(layer==11 && row==107) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==11 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=24; ++ch) if(layer==11 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=28; ch<=44; ++ch) if(layer==11 && row==52+(ch-28)*2) flag=1;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==11 && row==86+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==11 && row==106+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==10 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==10 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==10 && row==60+(ch-32)*2) flag=2;
    for(Int_t ch=49; ch<=54; ++ch) if(layer==10 && row==92+(ch-49)*2) flag=2;
    if(layer==10 && row==86) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==10 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=25; ++ch) if(layer==10 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=27; ch<=29; ++ch) if(layer==10 && row==51+(ch-27)*2) flag=3;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==10 && row==61+(ch-32)*2) flag=3;
    for(Int_t ch=48; ch<=54; ++ch) if(layer==10 && row==91+(ch-48)*2) flag=3;
    if(layer==10 && row==87) flag=3;
    return flag;
  case 11:
    //AGET-0
    for(Int_t ch=0; ch<=1; ++ch) if(layer==11 && row==110+ch*2) flag=0;
    for(Int_t ch=4; ch<=10; ++ch) if(layer==11 && row==118+(ch-4)*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==132+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==11 && row==152+(ch-23)*2) flag=0;
    for(Int_t ch=33; ch<=44; ++ch) if(layer==11 && row==172+(ch-33)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==196+(ch-46)*2) flag=0;
    if(layer==11 && row==216) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=1; ++ch) if(layer==11 && row==109+ch*2) flag=1;
    for(Int_t ch=5; ch<=10; ++ch) if(layer==11 && row==119+(ch-5)*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==11 && row==131+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=30; ++ch) if(layer==11 && row==151+(ch-23)*2) flag=1;
    for(Int_t ch=34; ch<=44; ++ch) if(layer==11 && row==173+(ch-34)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==11 && row==195+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==11 && row==215+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=5; ++ch) if(layer==10 && row==105+ch*2) flag=2;
    for(Int_t ch=8; ch<=10; ++ch) if(layer==10 && row==121+(ch-8)*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==127+(ch-12)*2) flag=2;
    for(Int_t ch=25; ch<=44; ++ch) if(layer==10 && row==151+(ch-25)*2) flag=2;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==10 && row==191+(ch-46)*2) flag=2;
    if(layer==10 && row==147) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=6; ++ch) if(layer==10 && row==104+ch*2) flag=3;
    for(Int_t ch=8; ch<=10; ++ch) if(layer==10 && row==120+(ch-8)*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==10 && row==126+(ch-12)*2) flag=3;
    for(Int_t ch=26; ch<=28; ++ch) if(layer==10 && row==152+(ch-26)*2) flag=3;
    for(Int_t ch=30; ch<=44; ++ch) if(layer==10 && row==160+(ch-30)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==10 && row==190+(ch-46)*2) flag=3;
    if(layer==10 && row==146) flag=3;
    return flag;
  case 12:
    //AGET-0
    for(Int_t ch=0; ch<=9; ++ch) if(layer==13 && row==1+ch*2) flag=0;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==13 && row==25+(ch-13)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==13 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==13 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=9; ++ch) if(layer==13 && row==ch*2) flag=1;
    for(Int_t ch=13; ch<=21; ++ch) if(layer==13 && row==24+(ch-13)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==13 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==13 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==12 && row==1+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==43+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==87+(ch-46)*2) flag=2;
    for(Int_t ch=58; ch<=59; ++ch) if(layer==12 && row==109+(ch-58)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==12 && row==ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==42+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==86+(ch-46)*2) flag=3;
    for(Int_t ch=58; ch<=59; ++ch) if(layer==12 && row==108+(ch-58)*2) flag=3;
    return flag;
  case 13:
    //AGET-0
    for(Int_t ch=1; ch<=10; ++ch) if(layer==13 && row==110+(ch-1)*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==13 && row==130+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=42; ++ch) if(layer==13 && row==150+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==13 && row==194+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=2; ch<=10; ++ch) if(layer==13 && row==111+(ch-2)*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==13 && row==129+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=43; ++ch) if(layer==13 && row==149+(ch-23)*2) flag=1;
    for(Int_t ch=47; ch<=55; ++ch) if(layer==13 && row==195+(ch-47)*2) flag=1;
    if(layer==13 && row==213) flag=1;
    //AGET-2
    for(Int_t ch=4; ch<=10; ++ch) if(layer==12 && row==124+(ch-4)*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==138+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==158+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==202+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=60; ++ch) if(layer==12 && row==222+(ch-57)*2) flag=2;
    if(layer==12 && row==118) flag=2;
    if(layer==12 && row==120) flag=2;
    //AGET-3
    for(Int_t ch=5; ch<=10; ++ch) if(layer==12 && row==125+(ch-5)*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==12 && row==137+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==12 && row==157+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==12 && row==201+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==12 && row==221+(ch-57)*2) flag=3;
    if(layer==12 && row==119) flag=3;
    if(layer==12 && row==121) flag=3;
    return flag;
  case 14:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==15 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==15 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==86+(ch-46)*2) flag=1;
    if(layer==15 && row==106) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==14 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==14 && row==86+(ch-46)*2) flag=2;
    for(Int_t ch=54; ch<=55; ++ch) if(layer==14 && row==102+(ch-54)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==14 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=50; ++ch) if(layer==14 && row==87+(ch-46)*2) flag=3;
    for(Int_t ch=54; ch<=55; ++ch) if(layer==14 && row==103+(ch-54)*2) flag=3;
    return flag;
  case 15:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==15 && row==108+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==130+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==150+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==194+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==15 && row==107+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==15 && row==129+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==15 && row==149+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==15 && row==193+(ch-46)*2) flag=1;
    if(layer==15 && row==213) flag=1;
    //AGET-2
    if(layer==14 && row==107) flag=2;
    if(layer==14 && row==109) flag=2;
    for(Int_t ch=4; ch<=10; ++ch) if(layer==14 && row==115+(ch-4)*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==129+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==149+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==14 && row==193+(ch-46)*2) flag=2;
    //AGET-3
    if(layer==14 && row==106) flag=3;
    if(layer==14 && row==108) flag=3;
    for(Int_t ch=5; ch<=10; ++ch) if(layer==14 && row==116+(ch-5)*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==14 && row==128+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==14 && row==148+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==14 && row==192+(ch-46)*2) flag=3;
    return flag;
  case 16:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==17 && row==ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==22+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==42+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==86+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==17 && row==1+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==23+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==43+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==87+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==16 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=49; ++ch) if(layer==16 && row==86+(ch-46)*2) flag=2;
    for(Int_t ch=52; ch<=55; ++ch) if(layer==16 && row==98+(ch-52)*2) flag=2;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==16 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==16 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=48; ++ch) if(layer==16 && row==87+(ch-46)*2) flag=3;
    for(Int_t ch=52; ch<=55; ++ch) if(layer==16 && row==99+(ch-52)*2) flag=3;
    if(layer==16 && row==10) flag=3;
    if(layer==16 && row==109) flag=3;
    return flag;
  case 17:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==17 && row==113+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==135+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==155+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==199+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==219+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==17 && row==112+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==17 && row==134+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==17 && row==154+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==17 && row==198+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=59; ++ch) if(layer==17 && row==218+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=5; ++ch) if(layer==16 && row==111+ch*2) flag=2;
    for(Int_t ch=8; ch<=10; ++ch) if(layer==16 && row==127+(ch-8)*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==133+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==153+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==16 && row==197+(ch-46)*2) flag=2;
    if(layer==16 && row==217) flag=2;
    if(layer==16 && row==219) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=5; ++ch) if(layer==16 && row==110+ch*2) flag=3;
    for(Int_t ch=9; ch<=10; ++ch) if(layer==16 && row==128+(ch-9)*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==16 && row==132+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==16 && row==152+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==16 && row==196+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==16 && row==216+(ch-57)*2) flag=3;
    return flag;
  case 18:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==19 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==87+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=62; ++ch) if(layer==19 && row==107+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==19 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==86+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==19 && row==106+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==18 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==86+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==18 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==87+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==107+(ch-57)*2) flag=3;
    return flag;
  case 19:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==19 && row==120+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==142+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==162+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==206+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=62; ++ch) if(layer==19 && row==226+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==19 && row==119+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==19 && row==141+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==19 && row==161+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==19 && row==205+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=63; ++ch) if(layer==19 && row==225+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==18 && row==117+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==139+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==159+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==203+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==223+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==18 && row==116+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==18 && row==138+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==18 && row==158+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==18 && row==202+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==18 && row==222+(ch-57)*2) flag=3;
    return flag;
  case 20:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==21 && row==ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==22+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==42+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==86+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==106+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==21 && row==1+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==23+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==43+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==87+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==107+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==20 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==86+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==106+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==20 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==87+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==107+(ch-57)*2) flag=3;
    return flag;
  case 21:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==21 && row==117+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==139+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==159+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==203+(ch-46)*2) flag=0;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==223+(ch-57)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==21 && row==116+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==21 && row==138+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==21 && row==158+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==21 && row==202+(ch-46)*2) flag=1;
    for(Int_t ch=57; ch<=61; ++ch) if(layer==21 && row==222+(ch-57)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==20 && row==123+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==145+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==165+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==209+(ch-46)*2) flag=2;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==229+(ch-57)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==20 && row==122+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==20 && row==144+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==20 && row==164+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==20 && row==208+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=64; ++ch) if(layer==20 && row==228+(ch-57)*2) flag=3;
    return flag;
  case 22:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==23 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==23 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==23 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==23 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==22 && row==1+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==22 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=35; ++ch) if(layer==22 && row==43+(ch-23)*2) flag=2;
    for(Int_t ch=39; ch<=44; ++ch) if(layer==22 && row==75+(ch-39)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==87+(ch-46)*2) flag=2;
    if(layer==22 && row==107) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==22 && row==ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==22 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=36; ++ch) if(layer==22 && row==42+(ch-23)*2) flag=3;
    for(Int_t ch=39; ch<=44; ++ch) if(layer==22 && row==74+(ch-39)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==86+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==22 && row==106+(ch-57)*2) flag=3;
    return flag;
  case 23:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==23 && row==106+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==128+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==148+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==23 && row==192+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==23 && row==105+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==23 && row==127+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==23 && row==147+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==23 && row==191+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==22 && row==110+ch*2) flag=2;
    for(Int_t ch=12; ch<=17; ++ch) if(layer==22 && row==132+(ch-12)*2) flag=2;
    if(layer==22 && row==150) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==22 && row==152+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==196+(ch-46)*2) flag=2;
    if(layer==22 && row==216) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==22 && row==109+ch*2) flag=3;
    for(Int_t ch=12; ch<=18; ++ch) if(layer==22 && row==131+(ch-12)*2) flag=3;
    if(layer==22 && row==149) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==22 && row==151+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==22 && row==195+(ch-46)*2) flag=3;
    for(Int_t ch=57; ch<=58; ++ch) if(layer==22 && row==215+(ch-57)*2) flag=3;
    return flag;
  case 24:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==25 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==43+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==25 && row==87+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==25 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==42+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==25 && row==86+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==24 && row==1+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==24 && row==43+(ch-23)*2) flag=2;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==24 && row==61+(ch-32)*2) flag=2;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==24 && row==87+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==24 && row==ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=29; ++ch) if(layer==24 && row==42+(ch-23)*2) flag=3;
    for(Int_t ch=33; ch<=44; ++ch) if(layer==24 && row==62+(ch-33)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==24 && row==86+(ch-46)*2) flag=3;
    return flag;
  case 25:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==25 && row==102+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==124+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==144+(ch-23)*2) flag=0;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==25 && row==188+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==25 && row==101+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==25 && row==123+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==25 && row==143+(ch-23)*2) flag=1;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==25 && row==187+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==24 && row==104+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==126+(ch-12)*2) flag=2;
    for(Int_t ch=26; ch<=44; ++ch) if(layer==24 && row==152+(ch-26)*2) flag=2;
    for(Int_t ch=46; ch<=53; ++ch) if(layer==24 && row==190+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==24 && row==103+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==24 && row==125+(ch-12)*2) flag=3;
    for(Int_t ch=26; ch<=44; ++ch) if(layer==24 && row==151+(ch-26)*2) flag=3;
    for(Int_t ch=46; ch<=54; ++ch) if(layer==24 && row==189+(ch-46)*2) flag=3;
    return flag;
  case 26:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==27 && row==ch*2) flag=0;
    for(Int_t ch=12; ch<=20; ++ch) if(layer==27 && row==22+(ch-12)*2) flag=0;
    for(Int_t ch=24; ch<=44; ++ch) if(layer==27 && row==44+(ch-24)*2) flag=0;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==86+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==27 && row==1+ch*2) flag=1;
    for(Int_t ch=12; ch<=19; ++ch) if(layer==27 && row==23+(ch-12)*2) flag=1;
    if(layer==27 && row==41) flag=1;
    for(Int_t ch=24; ch<=44; ++ch) if(layer==27 && row==45+(ch-24)*2) flag=1;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==87+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==26 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==86+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==26 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==87+(ch-46)*2) flag=3;
    return flag;
  case 27:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==27 && row==99+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==27 && row==121+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==27 && row==141+(ch-23)*2) flag=0;
    for(Int_t ch=31; ch<=44; ++ch) if(layer==27 && row==157+(ch-31)*2) flag=0;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==185+(ch-46)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==27 && row==98+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==27 && row==120+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=28; ++ch) if(layer==27 && row==140+(ch-23)*2) flag=1;
    if(layer==27 && row==154) flag=1;
    for(Int_t ch=32; ch<=44; ++ch) if(layer==27 && row==158+(ch-32)*2) flag=1;
    for(Int_t ch=46; ch<=51; ++ch) if(layer==27 && row==184+(ch-46)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==26 && row==101+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==123+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==143+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==187+(ch-46)*2) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==26 && row==100+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==26 && row==122+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==26 && row==142+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=52; ++ch) if(layer==26 && row==186+(ch-46)*2) flag=3;
    return flag;
  case 28:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==29 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==29 && row==43+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==29 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=34; ++ch) if(layer==29 && row==42+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==28 && row==1+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==23+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==43+(ch-23)*2) flag=2;
    if(layer==28 && row==87) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==28 && row==ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==22+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==42+(ch-23)*2) flag=3;
    if(layer==28 && row==86) flag=3;
    if(layer==28 && row==88) flag=3;
    return flag;
  case 29:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==29 && row==66+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==88+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=33; ++ch) if(layer==29 && row==108+(ch-23)*2) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==29 && row==65+ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==29 && row==87+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=34; ++ch) if(layer==29 && row==107+(ch-23)*2) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==28 && row==90+ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==112+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==132+(ch-23)*2) flag=2;
    if(layer==28 && row==176) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==28 && row==89+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==28 && row==111+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==28 && row==131+(ch-23)*2) flag=3;
    if(layer==28 && row==175) flag=3;
    if(layer==28 && row==177) flag=3;
    return flag;
  case 30:
    //AGET-0
    for(Int_t ch=0; ch<=10; ++ch) if(layer==31 && row==1+ch*2) flag=0;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==31 && row==23+(ch-12)*2) flag=0;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==31 && row==43+(ch-23)*2) flag=0;
    if(layer==31 && row==87) flag=0;
    if(layer==31 && row==89) flag=0;
    //AGET-1
    for(Int_t ch=0; ch<=10; ++ch) if(layer==31 && row==ch*2) flag=1;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==31 && row==22+(ch-12)*2) flag=1;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==31 && row==42+(ch-23)*2) flag=1;
    if(layer==31 && row==86) flag=1;
    if(layer==31 && row==88) flag=1;
    //AGET-2
    for(Int_t ch=0; ch<=10; ++ch) if(layer==30 && row==ch*2) flag=2;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==30 && row==22+(ch-12)*2) flag=2;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==30 && row==42+(ch-23)*2) flag=2;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==30 && row==86+(ch-46)*2) flag=2;
    if(layer==30 && row==106) flag=2;
    //AGET-3
    for(Int_t ch=0; ch<=10; ++ch) if(layer==30 && row==1+ch*2) flag=3;
    for(Int_t ch=12; ch<=21; ++ch) if(layer==30 && row==23+(ch-12)*2) flag=3;
    for(Int_t ch=23; ch<=44; ++ch) if(layer==30 && row==43+(ch-23)*2) flag=3;
    for(Int_t ch=46; ch<=55; ++ch) if(layer==30 && row==87+(ch-46)*2) flag=3;
    if(layer==30 && row==107) flag=3;
    return flag;
  default:
    return -1;
  }
}

//_____________________________________________________________________________
inline Int_t GetASADId(Int_t layer, Int_t row) //0~30
{
  Int_t flag=layer/4;
  Int_t section;
  if(flag==0) section=0; //layer 0~3
  else if(flag==1) section=1; //layer 4~7
  else if(layer==30||layer==31) section=3; //layer 30~31
  else section=2; //layer 8~29

  Int_t half=padParameter[layer][1]/2;
  Int_t division1=padParameter[layer][1]/6;
  Int_t division2=padParameter[layer][1]*5/6;

  switch(section){
  case 0:
    if(row<half)
      return 0;
    else
      return 1;
  case 1:
    Int_t dummy;
    if(layer%4<2) dummy=2;
    if(2<=layer%4) dummy=5;
    if(row<division1||division2<=row)
      return dummy;
    else if(division1<=row&&row<half)
      return dummy+1;
    else
      return dummy+2;
  case 3:
    return 30;
  default:
    if(row<half)
      return layer-layer%2;
    else
      return layer-layer%2+1;
  }
}

//_____________________________________________________________________________
inline Int_t GetCoBoId(Int_t layer, Int_t row)
{
  switch(layer){
  case 4:
    if(60<=row && row<=99)
      return 1;
    else
      return 0;
  case 5:
    if(72<=row && row<=119)
      return 1;
    else
      return 0;
  default:
    return layer/4;
  }
}

//_____________________________________________________________________________
inline Int_t GetPadId(Int_t layerID, Int_t rowID)
{
  //Original
  //Int_t padID=0;
  // Check!!!!!!!
  Int_t padID=1;
  for(Int_t layi = 0 ; layi<layerID; layi++) padID += padParameter[layi][1];
  padID+=rowID;
  return padID;

}

inline Int_t getLayerID(Int_t padID)
{
  padID-=1;
  Int_t layer;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  return layer;
}

inline Int_t getRowID(Int_t padID)
{
  padID-=1;
  Int_t layer, row;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;
  return row;
}
/*
  Double_t getTheta(Int_t layerID, Int_t rowID)
  {
  Double_t sTheta = 180.-(360./padParameter[layerID][3])*padParameter[layerID][1]/2.;
  Double_t theta = sTheta+(rowID+0.5)*(360.-2*sTheta)/padParameter[layerID][1];
  return theta;
  }
*/

inline Double_t getTheta(Int_t padID)
{
  padID-=1;
  Int_t layer, row;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;
  //std::cout<<"layer="<<layer<<", row="<<row<<std::endl;
  Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
  //Double_t theta = sTheta+(row+0.5)*(360.-2*sTheta)/padParameter[layer][1];
  Double_t theta = sTheta+(row+0.5)*360./padParameter[layer][3]-180;
  //std::cout<<"theta="<<theta<<std::endl;

  return theta;
}

inline Double_t getTheta(Int_t layer, Double_t m_row)
{
  Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
  //Double_t theta = sTheta+(row+0.5)*(360.-2*sTheta)/padParameter[layer][1];
  Double_t theta = sTheta+(m_row+0.5)*360./padParameter[layer][3]-180;
  //std::cout<<"theta="<<theta<<std::endl;
  return theta;
}

inline Double_t getMrow(Int_t layer, Double_t m_phi)
{

  Double_t mrow = 0.5*(padParameter[layer][1]-1.) + (90.-m_phi)*padParameter[layer][3]/360.;
  if(mrow<-0.0001){
    mrow = 0.5*(padParameter[layer][1]-1.) + (450.-m_phi)*padParameter[layer][3]/360.;
  }
  return mrow;
}

inline Double_t GetRadius(Int_t layer)
{
  return padParameter[layer][2];
}

inline Double_t getR(Int_t padID)
{
  padID-=1;
  Int_t layer;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  Double_t R = padParameter[layer][2];
  return R;
}

inline TVector3 getPosition(Int_t padID)
{
  padID-=1;
  Int_t layer, row;
  Int_t sum = 0;

  for (layer = 0; layer <= 30 && sum + padParameter[layer][1] <= padID; layer++)
  {
    sum += padParameter[layer][1];
  }
  row = padID - sum;

  TVector3 result;
  if (row > padParameter[layer][1]){ // out of range
    result.SetX(0);
    result.SetY(-1);
    result.SetZ(0);
  }
  else{
    Double_t x, z;
    //Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;

    //    x = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
    //    z = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) + ZTarget;

    // x = padParameter[layer][2] * sin(getTheta(padID+1)*TMath::Pi()/180.);
    // z = padParameter[layer][2] * cos(getTheta(padID+1)*TMath::Pi()/180.) + ZTarget;
    //std::cout<<"layer="<<layer<<", row"<<row<<std::endl;
    x = padParameter[layer][2] * sin(getTheta(layer,row)*TMath::Pi()/180.);
    z = padParameter[layer][2] * cos(getTheta(layer,row)*TMath::Pi()/180.) + ZTarget;

    // Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;
    // Double_t x_ = padParameter[layer][2] * -sin((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.);
    // Double_t z_ = padParameter[layer][2] * -cos((360./padParameter[layer][3])*TMath::Pi()/180. * (row + 0.5) + sTheta*TMath::Pi()/180.) + ZTarget;
    //std::cout<<"x="<<x<<", z="<<z<<", x_="<<x_<<", z_="<<z_<<std::endl;
    result.SetX(x);
    result.SetY(0);
    result.SetZ(z);
  }
  return result;
}

inline TVector3 getPosition(Int_t layer, Double_t m_row)
{
  TVector3 result;
  if(m_row > padParameter[layer][1]){ // out of range
    return TVector3(TMath::QuietNaN(), TMath::QuietNaN(), TMath::QuietNaN());
  }
  else{
    return TVector3(
      padParameter[layer][2] * sin(getTheta(layer, m_row)*TMath::Pi()/180.),
      0.,
      padParameter[layer][2] * cos(getTheta(layer, m_row)*TMath::Pi()/180.) + ZTarget);
  }
}

inline Int_t findPadID(Double_t z, Double_t x)
{
  z -= ZTarget;
  Double_t radius = sqrt(x*x + z*z);
  Double_t angle;
  if (z == 0)
  {
    if (x > 0)   angle = 1.5*TMath::Pi();
    else if (x < 0)   angle = 0.5*TMath::Pi();
    else return -1000; // no padID if (0,0)
  }
  //  else
  //  {
  //    if (z < 0 && x < 0) angle = atan(x / z);
  //    else if (z > 0 && x < 0) angle = TMath::Pi() - atan(-x / z);
  //    else if (z > 0 && x > 0) angle = TMath::Pi() + atan(x / z);
  //    else if (z < 0 && x > 0) angle = 2*TMath::Pi() - atan(-x / z);
  //  }
//  else if (z < 0) angle = atan(x / z);
//  else  angle = TMath::Pi() + atan(x / z);
  else{
		if (z > 0) angle = TMath::Pi()+atan(x / z);
	  else if( z < 0&&x<0) angle = atan(x / z);
		  else angle = 2*TMath::Pi()+ atan(x / z);//angle of z<0&&x>0 plane should be [1.5Pi,2Pi], not [-0.5Pi , 0].
	}
	//cout << " angle: " << angle*180/TMath::Pi() << endl;

  Int_t layer, row;
  // find layer_num.
  for (layer = 0; layer<NumOfLayersTPC;layer++)
  {
    if (layer != 0)
    {
      if (padParameter[layer][2] - padParameter[layer][5] * 0.5 >= radius &&
          padParameter[layer - 1][2] + padParameter[layer - 1][5] * 0.5 <= radius) return -layer;
    }
		double rad_in= padParameter[layer][2]-padParameter[layer][5]*0.5;
		double rad_out= padParameter[layer][2]+padParameter[layer][5]*0.5;
  	if(rad_in<=radius and rad_out>=radius){
			break;
		}
		if(layer==NumOfLayersTPC-1 && rad_out<radius) return -1000;
	}
  //cout << " layer: " << layer << endl;

  Double_t sTheta = 180.-(360./padParameter[layer][3])*padParameter[layer][1]/2.;

  // find row_num
  if (angle - (sTheta*TMath::Pi()/180.) < 0) return -1000;

  //Double_t a, b, c;
  row = (int)((angle-(sTheta*TMath::Pi()/180.))/(360./padParameter[layer][3]*TMath::Pi()/180.));
  if (row > padParameter[layer][1]) return -1000;

  //cout << " row: " << row << endl;
  //This is original one
  //return GetPadId(layer, row)+1;
  //Please check
  return GetPadId(layer, row);
}

//_____________________________________________________________________________
inline Double_t
ArcLength(Int_t layer, Double_t row1, Double_t row2)
{
  const Int_t R = padParameter[layer][2];
  Double_t theta = getTheta(layer, row1) - getTheta(layer, row2);
  theta = std::fmod(theta, 2*TMath::Pi());
  if(theta < 0) theta += 2*TMath::Pi();
  theta = TMath::Min(theta, 2*TMath::Pi() - theta);
  return R*theta;
}

//_____________________________________________________________________________
inline void
InitializeHistograms()
{
  std::vector<TH2Poly*> target;
  TList* list = gDirectory->GetList();
  TIter itr(list);
  while(itr.Next()){
    const TString& name((*itr)->GetName());
    const TString& cname((*itr)->ClassName());
    std::cout << " " << std::setw(8) << std::left << name
              << "(" << cname << ")" << std::endl;
    if(cname.EqualTo("TH2Poly")){
      target.push_back(dynamic_cast<TH2Poly*>(*itr));
    }
  }

  Double_t X[5];
  Double_t Y[5];

  for(Int_t layer=0; layer<NumOfLayersTPC; ++layer){
    Double_t pLength = tpc::padParameter[layer][5];
    Double_t st = 180.-(360./tpc::padParameter[layer][3])
      * tpc::padParameter[layer][1]/2.;
    Double_t sTheta  = (-1+st/180.)*TMath::Pi();
    Double_t dTheta  = (360./tpc::padParameter[layer][3])/180.*TMath::Pi();
    Double_t cRad    = tpc::padParameter[layer][2];
    Int_t    nPad    = tpc::padParameter[layer][1];
    for(Int_t j=0; j<nPad; ++j){
      X[1] = (cRad+(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[2] = (cRad+(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[3] = (cRad-(pLength/2.))*TMath::Cos((j+1)*dTheta+sTheta);
      X[4] = (cRad-(pLength/2.))*TMath::Cos(j*dTheta+sTheta);
      X[0] = X[4];
      Y[1] = (cRad+(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[2] = (cRad+(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[3] = (cRad-(pLength/2.))*TMath::Sin((j+1)*dTheta+sTheta);
      Y[4] = (cRad-(pLength/2.))*TMath::Sin(j*dTheta+sTheta);
      Y[0] = Y[4];
      for(Int_t ii=0; ii<5; ++ii) X[ii] += ZTarget;
      for(auto& h: target) h->AddBin(5, X, Y);
    }
  }
}

//_____________________________________________________________________________
inline Bool_t
IsClusterable(Int_t layer, Int_t row_a, Int_t row_b)
{

  Int_t deadpads=0;
  for(Int_t row=TMath::Min(row_a, row_b)+1;row<TMath::Max(row_a, row_b);row++){
    Int_t padID = GetPadId(layer, row);
    if(std::find(std::begin(deadChannel[layer]), std::end(deadChannel[layer]), padID) != std::end(deadChannel[layer])) deadpads++;
    if(std::find(std::begin(noisypad), std::end(noisypad), padID) != std::end(noisypad)) deadpads++;
  }

  if(layer < 10){
    const Int_t npad = padParameter[layer][kNumOfPad];
    return (TMath::Abs(row_a - row_b) <= MaxRowDifTPC + deadpads
            || TMath::Abs(row_a - row_b) >= npad - MaxRowDifTPC - deadpads);
  }else{
    return (TMath::Abs(row_a - row_b) <= MaxRowDifTPC + deadpads);
  }
}

//_____________________________________________________________________________
inline Bool_t Dead(Int_t padID){

  Bool_t centerframe = std::find(std::begin(padOnCenterFrame), std::end(padOnCenterFrame), padID) != std::end(padOnCenterFrame);
  Int_t layer = getLayerID(padID);
  Bool_t dead = std::find(std::begin(deadChannel[layer]), std::end(deadChannel[layer]), padID) != std::end(deadChannel[layer]);
  Bool_t noisy = std::find(std::begin(noisypad), std::end(noisypad), padID) != std::end(noisypad);
  if(centerframe||dead||noisy) return true;
  else return false;
}

//_____________________________________________________________________________
inline Bool_t Dead(Int_t layer, Int_t row){

  Int_t padID = GetPadId(layer, row);
  return Dead(padID);
}


//Functions for G4 simulation
inline Double_t GetClSize1Prob(double mom, int pid, int layer){
    bool inner = false;
    if(layer < 10) inner = true;
    int pidflag = 0;
    if(abs(pid)<300) pidflag = 0;
    else pidflag = 1;
    int mom_flag = (int)(mom/0.1);
    double prob = 0;
    if(mom_flag > 9){
      mom_flag = 9;
      if(inner){
        prob = ClusterSizeInner[pidflag][mom_flag];
      }
      else{
        prob = ClusterSizeOuter[pidflag][mom_flag];
      }
    }
    else{
      double resi = (mom )/0.1 - mom_flag;
      if(inner){
        prob = ClusterSizeInner[pidflag][mom_flag]*(1-resi) + ClusterSizeInner[pidflag][mom_flag+1]*resi; 
      }
      else{
        prob = ClusterSizeOuter[pidflag][mom_flag]*(1-resi) + ClusterSizeOuter[pidflag][mom_flag+1]*resi; 
      }
    }
    return prob; 
  }
  inline Double_t GetDetectionEfficiency(TVector3 pos, Int_t pid, TVector3 mom, Double_t de){
    Int_t pad = findPadID(pos.z(), pos.x());
    Int_t layer = getLayerID(pad);
    Int_t row = getRowID(pad);
    bool Inner = false;
    if(layer < 10)Inner = true;
    int pidflag = 0;// 0 -> pion, 1 -> kaon, 2 -> proton or any other baryons.
    if(abs(pid)== 211) pidflag = 0;
    else if(abs(pid)== 321) pidflag = 1;
    else if(abs(pid)== 2212 or abs(pid)> 3000) pidflag = 2;
    double eff = 0.;
    if(Inner){
      if(pidflag == 0 or pidflag == 1){
        eff = 1;
        if(layer>4){
          eff = 1 - 0.02*(layer-4);
        }
      }
      else if (pidflag == 2){
        eff = 1;
      }
    }
    else {
      if(pidflag == 0 or pidflag == 1){
        eff = 1;
        if(layer > 15 and layer < 29){
          eff = 1 - 0.1*(layer-15)/13;
        }
        else if (layer > 28 and layer < 31){
          eff = 0.9 - 0.05*(layer-28);
        }
        else if(layer > 30){
          eff = 0.8 + 0.1 *(layer - 30);
        }
      }
      else if (pidflag == 2){
        eff = 1;
      }
    }
    return eff;
  }
}
#endif

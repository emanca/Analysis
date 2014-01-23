#include "CMGTools/RootTools/interface/RochCor.h"
#include <TLorentzVector.h>

////^^^^^------------ GP BEGIN 
const double RochCor::pi = 3.14159265358979323846;
const float RochCor::genm_smr = 9.09956e+01; //gen mass peak with eta dependent gaussian smearing => better match in Z mass profile vs. eta/phi
const float RochCor::genm = 91.06; //gen mass peak without smearing => Z mass profile vs. eta/phi in CMS note
  
const float RochCor::recmA = 9.10062e+01; //rec mass peak in MC (2011A)
const float RochCor::drecmA = 9.09285e+01; //rec mass peak in data (2011A)
const float RochCor::mgsclA_stat = 0.0001; //stat. error of global factor for mass peak in MC (2011A)  
const float RochCor::mgsclA_syst = 0.0006; //syst. error of global factor for mass peak in MC (2011A)  
const float RochCor::dgsclA_stat = 0.0001; //stat. error of global factor for mass peak in data (2011A)
const float RochCor::dgsclA_syst = 0.0008; //syst. error of global factor for mass peak in data (2011A)
const float RochCor::recmB = 9.10210e+01; //rec mass peak in MC (2011B)
const float RochCor::drecmB = 9.09469e+01; //rec mass peak in data (2011B)
const float RochCor::mgsclB_stat = 0.0001; //stat. error of global factor for mass peak in MC (2011B)  
const float RochCor::mgsclB_syst = 0.0006; //syst. error of global factor for mass peak in MC (2011B)  
const float RochCor::dgsclB_stat = 0.0001; //stat. error of global factor for mass peak in data (2011B)
const float RochCor::dgsclB_syst = 0.0008; //syst. error of global factor for mass peak in data (2011B)
  
  //iteration2 after FSR : after Z Pt correction
const float RochCor::deltaA = -2.85242e-06;
const float RochCor::deltaA_stat = 7.74389e-07;
const float RochCor::deltaA_syst = 6.992e-07;
  
const float RochCor::sfA = 44.6463;
const float RochCor::sfA_stat = 1.92224;
const float RochCor::sfA_syst = 9.29;
  
const float RochCor::deltaB = -5.68463e-06;
const float RochCor::deltaB_stat = 8.21406e-07;
const float RochCor::deltaB_syst = 1.4268e-06;
  
const float RochCor::sfB = 23.8652;
const float RochCor::sfB_stat = 0.941748;
const float RochCor::sfB_syst = 4.86;

const float RochCor::apar = 1.0; //+- 0.002
const float RochCor::bpar = -5.03313e-06; //+- 1.57968e-06
const float RochCor::cpar = -4.41463e-05; //+- 1.92775e-06
const float RochCor::d0par = -0.000148871; //+- 3.16301e-06
const float RochCor::e0par = 1.59501; //+- 0.0249021
const float RochCor::d1par = 7.95495e-05; //+- 1.12386e-05
const float RochCor::e1par = -0.364823; //+- 0.17896
const float RochCor::d2par = 0.000152032; //+- 5.68386e-06
const float RochCor::e2par = 0.410195; //+- 0.0431732
////^^^^^------------ GP END 

const float RochCor::netabin[9] = {-2.4,-2.1,-1.4,-0.7,0.0,0.7,1.4,2.1,2.4};

const float RochCor::dcor_bfA[8][8]={{-0.000102967,-0.000025046,-0.000000182,-0.000031720,-0.000004638,-0.000013151,-0.000073829,-0.000021362},
				     {0.000075011,0.000054312,0.000003541,-0.000016074,-0.000013214,0.000000586,0.000025095,0.000117588},
				     {0.000147495,-0.000029203,0.000016442,0.000008401,-0.000014585,0.000004246,0.000027783,0.000023445},
				     {-0.000063689,-0.000021704,0.000006406,-0.000059618,-0.000025848,-0.000024249,-0.000044822,-0.000021290},
				     {-0.000000765,0.000011984,0.000027719,0.000025866,0.000017671,-0.000065838,-0.000047034,0.000044730},
				     {0.000011344,-0.000132266,-0.000038679,0.000015218,-0.000007268,-0.000022690,-0.000079248,-0.000052198},
				     {-0.000107277,-0.000092664,-0.000016977,-0.000022288,-0.000005622,-0.000042913,-0.000066225,-0.000058114},
				     {-0.000057816,-0.000107028,-0.000025582,0.000002045,-0.000035959,-0.000007281,-0.000059810,-0.000047769}};

const float RochCor::dcor_maA[8][8]={{0.000769889,0.000176340,-0.000173462,-0.000159710,-0.000081615,-0.000058009,0.000366711,0.001477802},
				     {0.001278711,0.000294977,-0.000105670,-0.000077729,-0.000068954,-0.000007808,0.000181101,0.000482557},
				     {0.000277706,0.000140310,0.000031660,0.000120208,0.000074286,0.000078156,0.000123767,0.000422373},
				     {-0.001571163,-0.000346042,0.000085722,0.000125968,0.000133283,0.000116788,0.000019394,0.000149867},
				     {-0.002099394,-0.000445580,0.000033833,0.000043528,0.000148554,0.000263179,0.000115391,-0.000513180},
				     {-0.001494752,-0.000433278,-0.000092362,-0.000026281,-0.000000523,0.000048183,0.000051742,-0.000317282},
				     {-0.000013164,-0.000104072,-0.000063807,-0.000056538,-0.000067794,-0.000125772,-0.000013945,0.001270347},
				     {0.000130353,0.000002891,-0.000136974,-0.000116878,-0.000190943,-0.000198251,0.000126934,0.001696517}};

const float RochCor::mcor_bfA[8][8]={{-0.000063713,-0.000029352,-0.000000867,0.000032270,0.000015492,0.000008083,-0.000069408,-0.000091716},
				     {-0.000060494,-0.000002986,0.000012797,-0.000031184,-0.000031340,-0.000006048,0.000013469,-0.000020202},
				     {-0.000022997,-0.000043807,-0.000007016,-0.000030670,-0.000020740,-0.000016735,-0.000007218,-0.000065682},
				     {-0.000041813,-0.000042280,0.000013533,-0.000002547,-0.000017769,-0.000011282,-0.000022693,-0.000099361},
				     {-0.000031962,-0.000022918,0.000009196,0.000027077,0.000002257,-0.000016681,-0.000017862,-0.000056932},
				     {-0.000026284,-0.000053526,-0.000000128,0.000026427,0.000034433,-0.000004638,-0.000023034,-0.000069140},
				     {-0.000109084,-0.000073483,-0.000007053,0.000037945,0.000037618,-0.000016044,-0.000053256,-0.000045541},
				     {-0.000063455,-0.000018084,-0.000009968,-0.000004891,-0.000018822,0.000001771,-0.000021826,-0.000079626}};

const float RochCor::mcor_maA[8][8]={{0.000950726,0.000132997,-0.000166230,-0.000178984,-0.000175606,-0.000184543,0.000028977,0.000145090},
				     {0.001082393,0.000012119,-0.000105033,-0.000095558,-0.000087842,-0.000050033,0.000203732,0.000781017},
				     {0.000522823,-0.000027809,0.000020088,0.000027120,0.000029425,0.000063659,0.000059290,0.000007311},
				     {0.000659471,0.000426387,0.000184802,0.000127485,0.000153550,0.000143188,0.000325332,0.000307829},
				     {0.000842162,0.000490264,0.000212897,0.000147332,0.000174670,0.000153595,0.000327076,0.001000893},
				     {-0.001242714,-0.000155280,0.000101135,0.000094522,0.000093880,0.000066729,0.000141144,0.000810823},
				     {-0.001757072,-0.000320008,-0.000029266,-0.000022502,-0.000040205,-0.000056041,-0.000149048,-0.000221401},
				     {0.000408788,0.000114598,-0.000141981,-0.000110819,-0.000115938,-0.000138071,-0.000784406,-0.002160131}};

const float RochCor::dcor_bfAer[8][8]={{0.000069681,0.000035377,0.000030160,0.000028300,0.000028481,0.000030308,0.000035908,0.000066403},
                                       {0.000063562,0.000035195,0.000029976,0.000028288,0.000028251,0.000030304,0.000035175,0.000063468},
                                       {0.000066084,0.000036266,0.000030191,0.000028046,0.000028118,0.000029895,0.000035353,0.000063882},
                                       {0.000064693,0.000035627,0.000029777,0.000028668,0.000028323,0.000030154,0.000034876,0.000063622},
                                       {0.000065655,0.000035484,0.000030380,0.000028062,0.000028263,0.000030324,0.000035823,0.000068903},
                                       {0.000062857,0.000034907,0.000029606,0.000028968,0.000028557,0.000029858,0.000034830,0.000063717},
                                       {0.000066211,0.000035707,0.000029803,0.000028436,0.000028707,0.000029851,0.000035014,0.000064730},
                                       {0.000065003,0.000035761,0.000030160,0.000028192,0.000028342,0.000029811,0.000035545,0.000063645}};

const float RochCor::dcor_maAer[8][8]={{0.000069681,0.000035377,0.000030160,0.000028300,0.000028481,0.000030308,0.000035908,0.000066403},
                                       {0.000063562,0.000035195,0.000029976,0.000028288,0.000028251,0.000030304,0.000035175,0.000063468},
                                       {0.000066084,0.000036266,0.000030191,0.000028046,0.000028118,0.000029895,0.000035353,0.000063882},
                                       {0.000064693,0.000035627,0.000029777,0.000028668,0.000028323,0.000030154,0.000034876,0.000063622},
                                       {0.000065655,0.000035484,0.000030380,0.000028062,0.000028263,0.000030324,0.000035823,0.000068903},
                                       {0.000062857,0.000034907,0.000029606,0.000028968,0.000028557,0.000029858,0.000034830,0.000063717},
                                       {0.000066211,0.000035707,0.000029803,0.000028436,0.000028707,0.000029851,0.000035014,0.000064730},
                                       {0.000065003,0.000035761,0.000030160,0.000028192,0.000028342,0.000029811,0.000035545,0.000063645}};

const float RochCor::mcor_bfAer[8][8]={{0.000028957,0.000015643,0.000013419,0.000012634,0.000012700,0.000013592,0.000016042,0.000028597},
				       {0.000027958,0.000015560,0.000013545,0.000012820,0.000012798,0.000013524,0.000015725,0.000027844},
				       {0.000027910,0.000015737,0.000013522,0.000012785,0.000012761,0.000013554,0.000015626,0.000027776},
				       {0.000028081,0.000015884,0.000013473,0.000012691,0.000012659,0.000013430,0.000015598,0.000027889},
				       {0.000027971,0.000015665,0.000013466,0.000012651,0.000012648,0.000013558,0.000016132,0.000029045},
				       {0.000027824,0.000015624,0.000013452,0.000012922,0.000012881,0.000013473,0.000015628,0.000027859},
				       {0.000028053,0.000015657,0.000013501,0.000012726,0.000012889,0.000013432,0.000015494,0.000027716},
				       {0.000028212,0.000015901,0.000013511,0.000012648,0.000012674,0.000013507,0.000015666,0.000027969}};

const float RochCor::mcor_maAer[8][8]={{0.000028957,0.000015643,0.000013419,0.000012634,0.000012700,0.000013592,0.000016042,0.000028597},
				      {0.000027958,0.000015560,0.000013545,0.000012820,0.000012798,0.000013524,0.000015725,0.000027844},
				      {0.000027910,0.000015737,0.000013522,0.000012785,0.000012761,0.000013554,0.000015626,0.000027776},
				      {0.000028081,0.000015884,0.000013473,0.000012691,0.000012659,0.000013430,0.000015598,0.000027889},
				      {0.000027971,0.000015665,0.000013466,0.000012651,0.000012648,0.000013558,0.000016132,0.000029045},
				      {0.000027824,0.000015624,0.000013452,0.000012922,0.000012881,0.000013473,0.000015628,0.000027859},
				      {0.000028053,0.000015657,0.000013501,0.000012726,0.000012889,0.000013432,0.000015494,0.000027716},
				      {0.000028212,0.000015901,0.000013511,0.000012648,0.000012674,0.000013507,0.000015666,0.000027969}};


//=======================================================================================================

const float RochCor::dmavgA[8][8]={{0.025922541,0.025094489,0.025024760,0.025459164,0.025507064,0.024926673,0.025264207,0.026154362},
				   {0.025771485,0.025052688,0.025031280,0.025448624,0.025418813,0.024947593,0.025213752,0.025711461},
				   {0.025992243,0.025246736,0.025081158,0.025465835,0.025396615,0.025090199,0.025225184,0.025674825},
				   {0.026065594,0.025210021,0.024985654,0.025468545,0.025506958,0.025011636,0.025137782,0.025900352},
				   {0.025723593,0.025225323,0.025057659,0.025327039,0.025445884,0.025039377,0.025250011,0.025920693},
				   {0.025890951,0.025184183,0.025108732,0.025431830,0.025389774,0.025015759,0.025133115,0.025839978},
				   {0.025969359,0.025120514,0.025090360,0.025397708,0.025439110,0.024991973,0.025145588,0.025956176},
				   {0.025890127,0.025152920,0.025013377,0.025419729,0.025451053,0.025052175,0.025151170,0.025819692}};

const float RochCor::dpavgA[8][8]={{0.025916064,0.025307681,0.024921634,0.025494383,0.025417695,0.025041578,0.025211546,0.026138035},
				   {0.025661559,0.025125726,0.024917529,0.025545957,0.025489507,0.024993156,0.025232820,0.025819186},
				   {0.025977549,0.025185495,0.024993264,0.025475434,0.025447391,0.025060103,0.025212464,0.025954195},
				   {0.026050003,0.025305680,0.024998627,0.025457135,0.025432215,0.025014203,0.025160101,0.025892301},
				   {0.026134513,0.025227871,0.024904078,0.025474157,0.025492296,0.025053230,0.025237275,0.025852681},
				   {0.026070419,0.025241236,0.025032270,0.025465782,0.025477130,0.025089989,0.025266698,0.025788486},
				   {0.025893934,0.025278195,0.025005032,0.025417408,0.025462482,0.025002584,0.025334630,0.025973266},
				   {0.025995267,0.025325988,0.024916915,0.025384798,0.025394341,0.025024826,0.025278556,0.025963080}};

const float RochCor::mmavgA[8][8]={{0.025859084,0.025171574,0.025031055,0.025414657,0.025421970,0.025005626,0.025202085,0.025905427},
				   {0.025820029,0.025133587,0.025017091,0.025474822,0.025465318,0.025012641,0.025166745,0.025890426},
				   {0.025843527,0.025164085,0.025038476,0.025478059,0.025477999,0.025038661,0.025147011,0.025803665},
				   {0.025876173,0.025204137,0.025021527,0.025462181,0.025430126,0.025017812,0.025204524,0.025839771},
				   {0.025861181,0.025205469,0.025018335,0.025426248,0.025460153,0.025022144,0.025196286,0.025888667},
				   {0.025758405,0.025144612,0.025037865,0.025474080,0.025474752,0.025036212,0.025167523,0.025862979},
				   {0.025807403,0.025158727,0.025014627,0.025464075,0.025472244,0.025021210,0.025130780,0.025783651},
				   {0.025905758,0.025201579,0.025009920,0.025423990,0.025416973,0.025032037,0.025147863,0.025760345}};

const float RochCor::mpavgA[8][8]={{0.025888818,0.025206935,0.025026146,0.025442136,0.025470197,0.025027800,0.025243407,0.025952872},
				   {0.025923013,0.025178692,0.025011691,0.025498907,0.025477111,0.025023833,0.025166403,0.025886905},
				   {0.025943006,0.025225756,0.025007080,0.025476595,0.025488864,0.025026658,0.025179579,0.025889859},
				   {0.025949377,0.025208820,0.025021340,0.025454777,0.025444410,0.024996581,0.025195810,0.025967101},
				   {0.025950823,0.025197201,0.024994630,0.025441154,0.025458403,0.025024723,0.025211150,0.025938079},
				   {0.025939334,0.025189405,0.025023127,0.025486532,0.025474880,0.025007298,0.025220661,0.025942054},
				   {0.025997086,0.025189140,0.025022336,0.025471484,0.025493603,0.025012869,0.025188219,0.025964187},
				   {0.025944809,0.025237280,0.025027926,0.025469723,0.025445395,0.025024964,0.025234456,0.026013442}};

//=======================================================================================================

const float RochCor::dcor_bfB[8][8]={{-0.000121996,-0.000051596,-0.000011541,-0.000024676,-0.000058614,0.000004092,-0.000042827,-0.000040838},
				     {-0.000082378,-0.000066345,-0.000047373,-0.000017058,-0.000021958,-0.000029804,-0.000053044,-0.000014822},
				     {0.000012351,-0.000031871,-0.000017504,-0.000017341,-0.000018939,-0.000036424,-0.000023220,-0.000006308},
				     {0.000019524,-0.000027703,-0.000010703,-0.000002277,-0.000061078,-0.000063380,-0.000063290,-0.000021621},
				     {-0.000052157,-0.000054030,-0.000039215,-0.000047173,-0.000021800,-0.000008816,-0.000041229,-0.000075721},
				     {-0.000059315,-0.000081392,-0.000015056,0.000009267,0.000015595,-0.000038434,-0.000008257,-0.000000816},
				     {-0.000011477,-0.000061045,-0.000023999,0.000018858,0.000002374,0.000010510,-0.000017883,0.000022914},
				     {-0.000070644,-0.000061816,-0.000037444,-0.000036912,0.000013680,-0.000003858,-0.000005998,-0.000005702}};

const float RochCor::dcor_maB[8][8]={{0.000648425,0.000184921,-0.000159733,-0.000209004,-0.000158411,-0.000118417,0.000343531,0.001274013},
				     {0.001107062,0.000215835,-0.000074671,-0.000078597,-0.000068598,-0.000038450,0.000119307,0.000441764},
				     {0.000299114,0.000047329,0.000028768,0.000041104,0.000068554,0.000057753,0.000062958,0.000434265},
				     {-0.001531225,-0.000369769,0.000034113,0.000106520,0.000164404,0.000166858,0.000026463,0.000183539},
				     {-0.001931603,-0.000423243,0.000014606,0.000077337,0.000176900,0.000297178,0.000124199,-0.000358460},
				     {-0.001242424,-0.000306289,-0.000036347,0.000069690,0.000008892,0.000084983,0.000127487,-0.000191898},
				     {0.000001016,-0.000111617,-0.000068977,-0.000044695,-0.000085620,-0.000113063,0.000016393,0.001409055},
				     {0.000351159,0.000019529,-0.000157996,-0.000187339,-0.000146448,-0.000186771,0.000173020,0.001701876}};

const float RochCor::mcor_bfB[8][8]={{-0.000072402,-0.000058879,0.000003018,0.000018626,0.000007212,-0.000005316,-0.000095954,-0.000125599},
				     {-0.000061846,-0.000003202,0.000014083,-0.000033514,-0.000039273,0.000004540,0.000009430,-0.000017489},
				     {-0.000031382,-0.000047138,-0.000019203,-0.000024139,-0.000038678,-0.000040859,0.000012155,-0.000065070},
				     {-0.000056882,-0.000031030,0.000024829,0.000013713,-0.000010394,-0.000020459,-0.000045276,-0.000057322},
				     {-0.000043675,0.000008459,0.000015752,0.000013816,-0.000008688,-0.000031616,-0.000032060,-0.000053715},
				     {0.000000487,-0.000064518,0.000019284,0.000045588,0.000028956,-0.000002070,-0.000029702,-0.000080368},
				     {-0.000108752,-0.000075013,-0.000016411,0.000050559,0.000042682,-0.000014198,-0.000064707,-0.000072518},
				     {-0.000039973,-0.000016213,0.000000336,0.000004890,-0.000014730,-0.000011327,-0.000006556,-0.000088452}};

const float RochCor:: mcor_maB[8][8]={{0.000974618,0.000128823,-0.000152265,-0.000186697,-0.000167048,-0.000204591,0.000008788,0.000152219},
				      {0.001104959,0.000007490,-0.000113518,-0.000115239,-0.000087118,-0.000054860,0.000210312,0.000761740},
				      {0.000495175,-0.000021281,0.000019935,0.000025729,0.000026502,0.000058051,0.000080846,0.000040006},
				      {0.000679721,0.000425605,0.000190407,0.000121199,0.000160504,0.000146688,0.000321539,0.000309533},
				      {0.000862263,0.000497438,0.000205395,0.000146674,0.000179833,0.000151782,0.000320139,0.000994401},
				      {-0.001264907,-0.000156277,0.000113008,0.000094040,0.000099709,0.000062566,0.000128300,0.000821240},
				      {-0.001737087,-0.000315566,-0.000028978,-0.000021989,-0.000026353,-0.000042322,-0.000148004,-0.000245085},
				      {0.000404257,0.000109031,-0.000120937,-0.000124930,-0.000110856,-0.000152148,-0.000770993,-0.002169186}};


const float RochCor::dcor_bfBer[8][8]={{0.000073392,0.000037888,0.000030941,0.000028542,0.000028942,0.000031034,0.000038166,0.000071498},
				       {0.000067919,0.000036855,0.000030957,0.000028460,0.000028904,0.000030892,0.000036525,0.000067771},
				       {0.000069890,0.000037252,0.000030693,0.000028233,0.000028295,0.000030924,0.000036761,0.000067732},
				       {0.000068865,0.000037515,0.000030440,0.000028521,0.000028758,0.000031061,0.000036545,0.000067957},
				       {0.000070287,0.000038039,0.000031099,0.000028541,0.000028535,0.000030889,0.000038192,0.000073368},
				       {0.000070529,0.000036023,0.000030708,0.000029111,0.000029156,0.000030485,0.000037271,0.000069426},
				       {0.000068987,0.000036834,0.000030454,0.000028355,0.000028894,0.000030568,0.000036321,0.000069434},
				       {0.000069238,0.000037352,0.000030916,0.000028682,0.000028282,0.000030943,0.000037054,0.000068690}};

const float RochCor::dcor_maBer[8][8]={{0.000073392,0.000037888,0.000030941,0.000028542,0.000028942,0.000031034,0.000038166,0.000071498},
				       {0.000067919,0.000036855,0.000030957,0.000028460,0.000028904,0.000030892,0.000036525,0.000067771},
				       {0.000069890,0.000037252,0.000030693,0.000028233,0.000028295,0.000030924,0.000036761,0.000067732},
				       {0.000068865,0.000037515,0.000030440,0.000028521,0.000028758,0.000031061,0.000036545,0.000067957},
				       {0.000070287,0.000038039,0.000031099,0.000028541,0.000028535,0.000030889,0.000038192,0.000073368},
				       {0.000070529,0.000036023,0.000030708,0.000029111,0.000029156,0.000030485,0.000037271,0.000069426},
				       {0.000068987,0.000036834,0.000030454,0.000028355,0.000028894,0.000030568,0.000036321,0.000069434},
				       {0.000069238,0.000037352,0.000030916,0.000028682,0.000028282,0.000030943,0.000037054,0.000068690}};

const float RochCor::mcor_bfBer[8][8]={{0.000031813,0.000016534,0.000013999,0.000013046,0.000013123,0.000014245,0.000016964,0.000031719},
				       {0.000030718,0.000016483,0.000014127,0.000013244,0.000013198,0.000014151,0.000016711,0.000030863},
				       {0.000030665,0.000016682,0.000014126,0.000013207,0.000013188,0.000014194,0.000016551,0.000030725},
				       {0.000030997,0.000016830,0.000014015,0.000013110,0.000013070,0.000014090,0.000016538,0.000030930},
				       {0.000030883,0.000016557,0.000014048,0.000013105,0.000013067,0.000014169,0.000017090,0.000032174},
				       {0.000030553,0.000016523,0.000014045,0.000013333,0.000013320,0.000014083,0.000016591,0.000030854},
				       {0.000030886,0.000016615,0.000014077,0.000013136,0.000013316,0.000014054,0.000016452,0.000030541},
				       {0.000031052,0.000016845,0.000014093,0.000013053,0.000013071,0.000014144,0.000016638,0.000031061}};

const float RochCor::mcor_maBer[8][8]={{0.000031813,0.000016534,0.000013999,0.000013046,0.000013123,0.000014245,0.000016964,0.000031719},
				       {0.000030718,0.000016483,0.000014127,0.000013244,0.000013198,0.000014151,0.000016711,0.000030863},
				       {0.000030665,0.000016682,0.000014126,0.000013207,0.000013188,0.000014194,0.000016551,0.000030725},
				       {0.000030997,0.000016830,0.000014015,0.000013110,0.000013070,0.000014090,0.000016538,0.000030930},
				       {0.000030883,0.000016557,0.000014048,0.000013105,0.000013067,0.000014169,0.000017090,0.000032174},
				       {0.000030553,0.000016523,0.000014045,0.000013333,0.000013320,0.000014083,0.000016591,0.000030854},
				       {0.000030886,0.000016615,0.000014077,0.000013136,0.000013316,0.000014054,0.000016452,0.000030541},
				       {0.000031052,0.000016845,0.000014093,0.000013053,0.000013071,0.000014144,0.000016638,0.000031061}};

//=======================================================================================================

const float RochCor::dmavgB[8][8]={{0.025938774,0.025266827,0.025022293,0.025374338,0.025372916,0.025005996,0.025240074,0.026067392},
				   {0.025894727,0.025191097,0.025013219,0.025419539,0.025390481,0.025019205,0.025198246,0.025935809},
				   {0.025912919,0.025224576,0.025030573,0.025419982,0.025415115,0.025043226,0.025213202,0.025937655},
				   {0.025903638,0.025244655,0.025023684,0.025383398,0.025399556,0.025025292,0.025254252,0.025945694},
				   {0.025855091,0.025245849,0.025018633,0.025392413,0.025386199,0.025032433,0.025258546,0.025938246},
				   {0.025898044,0.025210425,0.025039174,0.025411554,0.025418325,0.025034439,0.025230630,0.025873626},
				   {0.025908219,0.025191441,0.025016559,0.025382898,0.025416337,0.025011142,0.025181345,0.025906562},
				   {0.026026606,0.025235830,0.025005632,0.025369690,0.025375398,0.025027767,0.025254124,0.025961736}};

const float RochCor::dpavgB[8][8]={{0.025993525,0.025262074,0.025032936,0.025402927,0.025411082,0.025025205,0.025257010,0.025987816},
				   {0.026016457,0.025236788,0.025026049,0.025445870,0.025426455,0.025012605,0.025226581,0.026018325},
				   {0.026046038,0.025298641,0.025007364,0.025396512,0.025432880,0.025031159,0.025238138,0.025957892},
				   {0.026067569,0.025302059,0.025011169,0.025396687,0.025370463,0.024996758,0.025241644,0.026118577},
				   {0.026205140,0.025263103,0.025003300,0.025373572,0.025393818,0.025008397,0.025280542,0.026013350},
				   {0.026010447,0.025247521,0.025027158,0.025427902,0.025395642,0.025008772,0.025299819,0.026073297},
				   {0.026045245,0.025236978,0.025011242,0.025420223,0.025442011,0.025030245,0.025237170,0.025985402},
				   {0.026017886,0.025261768,0.025024263,0.025412894,0.025397047,0.025037108,0.025256271,0.025938284}};

const float RochCor::mmavgB[8][8]={{0.025922899,0.025253554,0.025024570,0.025376081,0.025366729,0.025000526,0.025238491,0.026010632},
				   {0.025886079,0.025183953,0.025007207,0.025415872,0.025388657,0.025019143,0.025195823,0.025930482},
				   {0.025934619,0.025227622,0.025035391,0.025418887,0.025417957,0.025042034,0.025207412,0.025915978},
				   {0.025975643,0.025271630,0.025027551,0.025387821,0.025396539,0.025023647,0.025264276,0.025925383},
				   {0.025951465,0.025261106,0.025024150,0.025395296,0.025387278,0.025029785,0.025264818,0.025972973},
				   {0.025864186,0.025203249,0.025038447,0.025414974,0.025419505,0.025033609,0.025228557,0.025950473},
				   {0.025857412,0.025190903,0.025020008,0.025387552,0.025414808,0.025008508,0.025180904,0.025854423},
				   {0.026021465,0.025234665,0.025004645,0.025369725,0.025376778,0.025027170,0.025226526,0.025825617}};

const float RochCor::mpavgB[8][8]={{0.025972655,0.025257853,0.025030870,0.025399462,0.025417248,0.025026830,0.025266203,0.026020858},
				   {0.026012686,0.025235199,0.025023719,0.025450310,0.025427728,0.025011407,0.025226418,0.026002827},
				   {0.026042418,0.025301180,0.025007263,0.025397603,0.025435786,0.025029248,0.025229272,0.025961233},
				   {0.025988419,0.025271262,0.025005581,0.025394773,0.025372960,0.024995658,0.025227015,0.026099192},
				   {0.026098501,0.025217605,0.024991889,0.025370425,0.025396433,0.025007798,0.025271868,0.025980209},
				   {0.026023494,0.025241938,0.025023119,0.025431178,0.025393405,0.025007651,0.025298716,0.026031118},
				   {0.026105678,0.025244646,0.025012836,0.025423016,0.025446063,0.025025461,0.025247461,0.026041311},
				   {0.026005132,0.025261163,0.025024802,0.025413312,0.025401136,0.025030499,0.025288682,0.026116802}};

//===============================================================================================
//parameters for Z pt correction

const float RochCor::ptlow[85] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5,
				  6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5,
				  10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5,
				  15.0, 15.5, 16.0, 16.5, 17.0, 17.5, 18.0, 18.5, 19.0, 19.5,
				  20.0, 20.5, 21.0, 21.5, 22.0, 22.5, 23.0, 23.5, 24.0, 24.5,
				  25.0, 26.0, 27.0, 28.0, 29.0,
				  30.0, 32.0, 34.0, 36.0, 38.0,
				  40.0, 44.0, 48.0, 52.0, 56.0,
				  60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0, 95.0,
				  100.0, 110.0, 120.0, 130.0, 140.0, 150.0, 175.0,
				  200.0, 250.0, 350.0, 500.0, 1000.0};

//int nptbins( sizeof(ptlow)/sizeof(float) - 1 );
  
const float RochCor::zptscl[84] = {1.49177,1.45654,1.36283,1.28569,1.2418,1.12336,1.10416,1.08731,0.994051,0.96532,
				   0.94427,0.932725,0.918082,0.899665,0.898398,0.927687,0.908047,0.892392,0.924027,0.945895,
				   0.937149,0.923983,0.923387,0.955362,0.947812,0.962943,0.948781,0.961555,0.95222,0.999207,
				   0.973884,0.993013,0.953487,0.951402,0.985583,0.986603,0.981388,1.00022,1.0294,0.964748,
				   0.974592,1.01546,0.992343,1.00101,0.990866,0.98982,1.02924,1.02265,0.967695,1.02411,
				   0.97331,1.01052,1.01561,0.992594,0.976504,1.01205,0.981111,1.00078,1.02078,1.00719,
				   1.0099,1.02865,1.03845,1.03254,1.09815,1.10263,1.06302,1.0725,1.14703,1.10574,
				   1.13911,1.16947,1.1709,1.11413,1.28793,1.18953,1.20212,1.18112,1.25471,1.15329,
				   1.14276,1.17223,1.09173,2.00229};
  
const float RochCor::zptscler[84] = {0.0270027,0.0154334,0.0115338,0.00958085,0.0084683,0.00736665,0.0069567,0.00671434,
				     0.00617693,0.00601943,0.00594735,0.00594569,0.00594903,0.00595495,0.00608115,0.00633704,
				     0.0063916,0.0064468,0.00678106,0.00706769,0.00717517,0.00727958,0.00747182,0.00785544,
				     0.00798754,0.00828787,0.00839147,0.00865826,0.00876775,0.00933276,0.00935768,0.0097289,
				     0.00962058,0.00983828,0.0103044,0.0104871,0.0106575,0.0110388,0.0114986,0.0111494,
				     0.0115202,0.0121059,0.0121345,0.0124923,0.0125972,0.0128401,0.0134519,0.0136279,
				     0.0133414,0.014186,0.00992195,0.0105984,0.0109484,0.0111756,0.0114579,0.00870013,
				     0.00904749,0.00970734,0.0104583,0.0109818,0.00837852,0.00939894,0.010415,0.0113433,
				     0.013007,0.0128788,0.0140174,0.0156993,0.0181717,0.019765,0.0222326,0.0249408,
				     0.0272806,0.0211706,0.0278087,0.0306654,0.0361387,0.041327,0.0341513,0.0440116,
				     0.0473006,0.0680212,0.149162,0.56279};

RochCor::~RochCor(){
}

RochCor::RochCor(){
  
  eran.SetSeed(123456);
  sran.SetSeed(1234);
  
  for(int i=0; i<8; ++i){
      for(int j=0; j<8; ++j){
          mptsys_mc_dm[i][j]=0;
          mptsys_mc_da[i][j]=0;
          mptsys_da_dm[i][j]=0;
          mptsys_da_da[i][j]=0;
      }
  }

}

RochCor::RochCor(int seed){
  eran.SetSeed(123456);
  sran.SetSeed(seed);

  for(int i=0; i<8; ++i){
      for(int j=0; j<8; ++j){
          mptsys_mc_dm[i][j]=sran.Gaus(0.0, 1.0);
          mptsys_mc_da[i][j]=sran.Gaus(0.0, 1.0);
          mptsys_da_dm[i][j]=sran.Gaus(0.0, 1.0);
          mptsys_da_da[i][j]=sran.Gaus(0.0, 1.0);
      }
  }
}

void RochCor::momcor_mc( TLorentzVector& mu, float charge, float sysdev, int runopt){
  
  //sysdev == num : deviation = num

  float ptmu = mu.Pt();
  float muphi = mu.Phi();
  float mueta = mu.Eta(); // same with mu.Eta() in Root

  float px = mu.Px();
  float py = mu.Py();
  float pz = mu.Pz();
  float e = mu.E();

  int mu_phibin = phibin(muphi);
  int mu_etabin = etabin(mueta);

  //float mptsys = sran.Gaus(0.0,sysdev);
  
  float dm = 0.0;
  float da = 0.0;
  
  if(runopt == 0){
    dm = (mcor_bfA[mu_phibin][mu_etabin] + mptsys_mc_dm[mu_phibin][mu_etabin]*mcor_bfAer[mu_phibin][mu_etabin])/mmavgA[mu_phibin][mu_etabin];
    da = mcor_maA[mu_phibin][mu_etabin] + mptsys_mc_da[mu_phibin][mu_etabin]*mcor_maAer[mu_phibin][mu_etabin];
  }else if(runopt == 1){
    dm = (mcor_bfB[mu_phibin][mu_etabin] + mptsys_mc_dm[mu_phibin][mu_etabin]*mcor_bfBer[mu_phibin][mu_etabin])/mmavgB[mu_phibin][mu_etabin];
    da = mcor_maB[mu_phibin][mu_etabin] + mptsys_mc_da[mu_phibin][mu_etabin]*mcor_maBer[mu_phibin][mu_etabin];
  }
  
  float cor = 1.0/(1.0 + dm + charge*da*ptmu);
  
  //for the momentum tuning - eta,phi,Q correction
  px *= cor;
  py *= cor;
  pz *= cor;
  e  *= cor;
  
  float recm = 0.0;
  //  float drecm = 0.0; 
  float delta = 0.0;
  float sf = 0.0;

  float gscler = 0.0;
  float deltaer = 0.0;
  float sfer = 0.0;
  
  if(runopt==0){
    recm = recmA;
    //   drecm = drecmA;
    
    delta = deltaA;
    sf = sfA;
    
    gscler = TMath::Sqrt( TMath::Power(mgsclA_stat,2) + TMath::Power(mgsclA_syst,2) );
    deltaer = TMath::Sqrt( TMath::Power(deltaA_stat,2) + TMath::Power(deltaA_syst,2) );
    sfer = TMath::Sqrt( TMath::Power(sfA_stat,2) + TMath::Power(sfA_syst,2) );
  }else if(runopt==1){
    recm = recmB;
    //    drecm = drecmB;
    
    delta = deltaB;
    sf = sfB;
    
    gscler = TMath::Sqrt( TMath::Power(mgsclB_stat,2) + TMath::Power(mgsclB_syst,2) );
    deltaer = TMath::Sqrt( TMath::Power(deltaB_stat,2) + TMath::Power(deltaB_syst,2) );
    sfer = TMath::Sqrt( TMath::Power(sfB_stat,2) + TMath::Power(sfB_syst,2) );
  }
  
  float tune = 1.0/(1.0 + (delta + sysdev*deltaer)*sqrt(px*px + py*py)*eran.Gaus(1.0,(sf + sysdev*sfer)));
  
  px *= (tune); 
  py *= (tune);  
  pz *= (tune);  
  e  *= (tune);   
      
  float gscl = (genm_smr/recm);
  
  px *= (gscl + sysdev*gscler);
  py *= (gscl + sysdev*gscler);
  pz *= (gscl + sysdev*gscler);
  e  *= (gscl + sysdev*gscler);
  
  mu.SetPxPyPzE(px,py,pz,e);
  
}


void RochCor::momcor_data( TLorentzVector& mu, float charge, float sysdev, int runopt){
  
  float ptmu = mu.Pt();

  float muphi = mu.Phi();
  float mueta = mu.Eta(); // same with mu.Eta() in Root

  float px = mu.Px();
  float py = mu.Py();
  float pz = mu.Pz();
  float e = mu.E();
  
  int mu_phibin = phibin(muphi);
  int mu_etabin = etabin(mueta);
  
  //float mptsys1 = sran.Gaus(0.0,sysdev);
  
  float dm = 0.0;
  float da = 0.0;
  
  if(runopt==0){
    dm = (dcor_bfA[mu_phibin][mu_etabin] + mptsys_da_dm[mu_phibin][mu_etabin]*dcor_bfAer[mu_phibin][mu_etabin])/dmavgA[mu_phibin][mu_etabin];
    da = dcor_maA[mu_phibin][mu_etabin] + mptsys_da_da[mu_phibin][mu_etabin]*dcor_maAer[mu_phibin][mu_etabin];
  }else if(runopt==1){
    dm = (dcor_bfB[mu_phibin][mu_etabin] + mptsys_da_dm[mu_phibin][mu_etabin]*dcor_bfBer[mu_phibin][mu_etabin])/dmavgB[mu_phibin][mu_etabin];
    da = dcor_maB[mu_phibin][mu_etabin] + mptsys_da_da[mu_phibin][mu_etabin]*dcor_maBer[mu_phibin][mu_etabin];
  }
  
  float cor = 1.0/(1.0 + dm + charge*da*ptmu);
  
  px *= cor;
  py *= cor;
  pz *= cor;
  e  *= cor;
  
  //after Z pt correction
  float recm = 0.0;
  float gscler = 0.0;
  
  if(runopt==0){
    recm = drecmA;
    gscler = TMath::Sqrt( TMath::Power(dgsclA_stat,2) + TMath::Power(dgsclA_syst,2) );
  }else if(runopt==1){
    recm = drecmB;
    gscler = TMath::Sqrt( TMath::Power(dgsclB_stat,2) + TMath::Power(dgsclB_syst,2) );
  }
  
  float gscl = (genm_smr/recm);
  
  px *= (gscl + sysdev*gscler);
  py *= (gscl + sysdev*gscler);
  pz *= (gscl + sysdev*gscler);
  e  *= (gscl + sysdev*gscler);
  
  mu.SetPxPyPzE(px,py,pz,e);
  
}

void RochCor::musclefit_data( TLorentzVector& mu, TLorentzVector& mubar){

  float dpar1 = 0.0;
  float dpar2 = 0.0;
  float epar1 = 0.0;
  float epar2 = 0.0;
  
  if(fabs(mu.PseudoRapidity())<=0.9){
    dpar1 = d0par;
    epar1 = e0par;
  }else if(mu.PseudoRapidity()>0.9){
    dpar1 = d1par;
    epar1 = e1par;
  }else if(mu.PseudoRapidity()<-0.9){
    dpar1 = d2par;
    epar1 = e2par;
  }

  if(fabs(mubar.PseudoRapidity())<=0.9){
    dpar2 = d0par;
    epar2 = e0par;
  }else if(mubar.PseudoRapidity()>0.9){
    dpar2 = d1par;
    epar2 = e1par;
  }else if(mubar.PseudoRapidity()<-0.9){
    dpar2 = d2par;
    epar2 = e2par;
  }

  float corr1 = 1.0 + bpar*mu.Pt() + (-1.0)*cpar*mu.Pt()*TMath::Sign(float(1.0),float(mu.PseudoRapidity()))*TMath::Power(mu.PseudoRapidity(),2)
    + (-1.0)*dpar1*mu.Pt()*sin(mu.Phi() + epar1);
  float corr2 = 1.0 + bpar*mubar.Pt() + (1.0)*cpar*mubar.Pt()*TMath::Sign(float(1.0),float(mubar.PseudoRapidity()))*TMath::Power(mubar.PseudoRapidity(),2)
    + (1.0)*dpar2*mubar.Pt()*sin(mubar.Phi() + epar2);
  
  float px1 = mu.Px();
  float py1 = mu.Py();
  float pz1 = mu.Pz();
  float e1 = mu.E();
  
  float px2 = mubar.Px();
  float py2 = mubar.Py();
  float pz2 = mubar.Pz();
  float e2 = mubar.E();

  px1 *= corr1;
  py1 *= corr1;
  pz1 *= corr1;
  e1 *= corr1;
  
  px2 *= corr2;
  py2 *= corr2;
  pz2 *= corr2;
  e2 *= corr2;
  
  mu.SetPxPyPzE(px1,py1,pz1,e1);
  mubar.SetPxPyPzE(px2,py2,pz2,e2);

}

Int_t RochCor::phibin(float phi){
  
  int nphibin = -1;
  
  for(int i=0; i<8; i++){
    if(-pi+(2.0*pi/8.0)*i <= phi && -pi+(2.0*pi/8.0)*(i+1) > phi){
      nphibin = i;
      break;
    }
  }
  
  return nphibin;
}

Int_t RochCor::etabin(float eta){

  int nbin = -1;
  
  for(int i=0; i<8; i++){
    if(netabin[i] <= eta && netabin[i+1] > eta){
      nbin = i;
      break;
    }
  }
  
  return nbin;
}

float RochCor::zptcor(float gzpt) {
  int ibin( 0 );
  
  // mcptscl[] = 84 bins: [0] and [83] are the underflow and overflow
  if ( gzpt > ptlow[nptbins] ) return nptbins-1;
  if ( gzpt < ptlow[0      ] ) return 0;
  
  for ( int i=0; i<nptbins; ++i ) {
    if ( gzpt>=ptlow[i] && gzpt<ptlow[i+1] ) { ibin=i; break; }
  }

  float zptwt = zptscl[ibin];

  return zptwt;
}

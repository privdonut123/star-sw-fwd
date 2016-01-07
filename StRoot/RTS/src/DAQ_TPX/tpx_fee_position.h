#ifndef _TPX_FEE_POSITION_

#define _TPX_FEE_POSITION_

/*
	New FEE position after the 2015 FEE remap
*/

static const unsigned char tpx_fee_position[6][36] = {

/* RDO location 1 */
	{162, 163, 169, 167, 166, 172, 170, 255, 255, 
	 255, 255, 255, 175, 180, 177, 159, 181, 160,
	 178, 179, 255, 255, 255, 255, 164, 165, 173, 
	 171, 168, 176, 174, 255, 255, 255, 255, 255,},

/* RDO location 2 */
	{138, 134, 154, 133, 145, 128, 143, 139, 158, 
	 127, 150, 129, 146, 135, 148, 155, 147, 156, 
	 151, 140, 152, 255, 255, 255, 157, 136, 142, 
	 130, 149, 132, 161, 141, 144, 131, 153, 137,},

/* RDO location 3 */
        {113, 100,  99, 114, 101, 120, 107, 106, 121,
         108, 255, 255, 115, 102, 117, 116, 122, 109,
         123, 255, 255, 255, 255, 255, 104, 119, 103,
         118, 105, 111, 125, 110, 124, 112, 255, 255}, 

/* RDO location 4 */
        { 85,  70,  69,  86,  71,  92,  78,  77, 93,
          79, 255, 255,  87,  89,  72,  88,  73, 94,
          96,  80,  95,  81, 255, 255,  75,  91, 74, 
          90,  76,  83,  98,  82,  97,  84, 255, 255},

/* RDO location 5 */
        { 53,  37,  55,  36,  54,  38,  61,  46,  63,
          45,  62,  47,  56,  40,  57,  39,  41,  64, 
          65,  48,  49, 255, 255, 255,  58,  43,  60, 
          42,  59,  44,  66,  51,  68,  50,  67,  52},

/* RDO location 6 */
        { 18,   1,  20,   0,  19,   2,  27,  10,  29,
           9,  28,  11,  21,   4,  23,   3,  22,   5,
          30,  13,  32,  12,  31,  14,  24,   7,  26,
           6,  25,   8,  33,  16,  35,  15,  34,  17}};


#endif
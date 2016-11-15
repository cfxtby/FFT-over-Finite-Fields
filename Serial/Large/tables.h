uint64_t *table64;
uint64_t *table32;
uint64_t *table16;
uint64_t *table8;
uint64_t *table4;
uint64_t *table2;

primitivePolynomial *p64;
primitivePolynomial *p32;
primitivePolynomial *p16;
primitivePolynomial *p8;
primitivePolynomial *p4;
primitivePolynomial *p2;

uint64_t pol64;
uint64_t pol32;
uint64_t pol16;
uint64_t pol8;
uint64_t pol4;
uint64_t pol2;

void buildPols() {

	p64 = getPrimitive();
	p32 = getPrimitive();
	p16 = getPrimitive();
	p8 = getPrimitive();
	p4 = getPrimitive();
	p2 = getPrimitive();

	pol64 = 0;
	pol32 = 0;
	pol16 = 0;
	pol8 = 0;
	pol4 = 0;
	pol2 = 0;

	setBit(&pol64,33);
	setBit(&pol64,30);
	setBit(&pol64,26);
	setBit(&pol64,25);
	setBit(&pol64,24);
	setBit(&pol64,23);
	setBit(&pol64,22);
	setBit(&pol64,21);
	setBit(&pol64,20);
	setBit(&pol64,18);
	setBit(&pol64,13);
	setBit(&pol64,12);
	setBit(&pol64,11);
	setBit(&pol64,10);
	setBit(&pol64,7);
	setBit(&pol64,5);
	setBit(&pol64,4);
	setBit(&pol64,2);
	setBit(&pol64,1);
	setBit(&pol64,0);

	insertPrimitiveElement(p64,64);
	insertPrimitiveElement(p64,33);
	insertPrimitiveElement(p64,30);
	insertPrimitiveElement(p64,26);
	insertPrimitiveElement(p64,25);
	insertPrimitiveElement(p64,24);
	insertPrimitiveElement(p64,23);
	insertPrimitiveElement(p64,22);
	insertPrimitiveElement(p64,21);
	insertPrimitiveElement(p64,20);
	insertPrimitiveElement(p64,18);
	insertPrimitiveElement(p64,13);
	insertPrimitiveElement(p64,12);
	insertPrimitiveElement(p64,11);
	insertPrimitiveElement(p64,10);
	insertPrimitiveElement(p64,7);
	insertPrimitiveElement(p64,5);
	insertPrimitiveElement(p64,4);
	insertPrimitiveElement(p64,2);
	insertPrimitiveElement(p64,1);
	insertPrimitiveElement(p64,0); 

	setBit(&pol32,15);
	setBit(&pol32,9);
	setBit(&pol32,7);
	setBit(&pol32,4);
	setBit(&pol32,3);
	setBit(&pol32,0);

	insertPrimitiveElement(p32,32);
	insertPrimitiveElement(p32,15);
	insertPrimitiveElement(p32,9);
	insertPrimitiveElement(p32,7);
	insertPrimitiveElement(p32,4);
	insertPrimitiveElement(p32,3);
	insertPrimitiveElement(p32,0);

	setBit(&pol16,5);
	setBit(&pol16,3);
	setBit(&pol16,2);
	setBit(&pol16,0);

	insertPrimitiveElement(p16,16);
	insertPrimitiveElement(p16,5);
	insertPrimitiveElement(p16,3);
	insertPrimitiveElement(p16,2);
	insertPrimitiveElement(p16,0);

	setBit(&pol8,4);
	setBit(&pol8,3);
	setBit(&pol8,2);
	setBit(&pol8,0);

	insertPrimitiveElement(p8,8);
	insertPrimitiveElement(p8,4);
	insertPrimitiveElement(p8,3);
	insertPrimitiveElement(p8,2);
	insertPrimitiveElement(p8,0);

	setBit(&pol4,1);
	setBit(&pol4,0);

	insertPrimitiveElement(p4,4);
	insertPrimitiveElement(p4,1);
	insertPrimitiveElement(p4,0);

	setBit(&pol2,1);
	setBit(&pol2,0);

	insertPrimitiveElement(p2,2);
	insertPrimitiveElement(p2,1);
	insertPrimitiveElement(p2,0);

}

primitivePolynomial* getPol(int m) {

	if(m==64) return p64;
	else if(m==32) return p32;
	else if(m==16) return p16;
	else if(m==8) return p8;
	else if(m==4) return p4;
	else return p2;

}

uint64_t getPrimPol(int m) {

	if(m==64) return pol64;
	else if(m==32) return pol32;
	else if(m==16) return pol16;
	else if(m==8) return pol8;
	else if(m==4) return pol4;
	else return pol2;

}

void freePols() {
	free(p64->terms);
	free(p64);
	free(p32->terms);
	free(p32);
	free(p16->terms);
	free(p16);
	free(p8->terms);
	free(p8);
	free(p4->terms);
	free(p4);
	free(p2->terms);
	free(p2);
}

void buildTable() {

table64 = (uint64_t *)malloc(sizeof(uint64_t)*64);

table64[63] = 15588707484923141439;
table64[62] = 8904600567183668890;
table64[61] = 1987634009112634340;
table64[60] = 16357799454061364322;
table64[59] = 15874185344704082069;
table64[58] = 9346551186767583484;
table64[57] = 12201470279193741501;
table64[56] = 2003542142868985263;
table64[55] = 15997637896822489273;
table64[54] = 13850347089535448476;
table64[53] = 1369242815936580601;
table64[52] = 18134003619941989084;
table64[51] = 6699905841332793031;
table64[50] = 1273743758326378506;
table64[49] = 10290523665025904390;
table64[48] = 17908228910333867766;
table64[47] = 15025579193700632411;
table64[46] = 7464063630536342018;
table64[45] = 15363497627494193475;
table64[44] = 3901388236763721102;
table64[43] = 6817160897129051684;
table64[42] = 18119765278828351118;
table64[41] = 17824810300707909138;
table64[40] = 7783638453117792155;
table64[39] = 5661579920225750636;
table64[38] = 4347884501264900129;
table64[37] = 10633681636318195913;
table64[36] = 12946368814566557862;
table64[35] = 6959859988055258630;
table64[34] = 8879768586127329464;
table64[33] = 7382463516588101015;
table64[32] = 4125863343146719424;
table64[31] = 6590804557738356910;
table64[30] = 7228336178356203901;
table64[29] = 13526989581464352498;
table64[28] = 13227140292500532871;
table64[27] = 17978339156933666767;
table64[26] = 2871183876035159548;
table64[25] = 15027438019396783516;
table64[24] = 5930438465407464505;
table64[23] = 12502290484167703138;
table64[22] = 8638890845374011113;
table64[21] = 4356953267644119759;
table64[20] = 11593005521033879483;
table64[19] = 11728183960584221855;
table64[18] = 11796744625775092113;
table64[17] = 337546862542194462;
table64[16] = 15079845503333857856;
table64[15] = 16964771062782307951;
table64[14] = 15361140383494460467;
table64[13] = 12581810914998629814;
table64[12] = 4191380083765225707;
table64[11] = 6484635769251602160;
table64[10] = 4506624024319379718;
table64[9] = 3695448190228955447;
table64[8] = 14887336961718343758;
table64[7] = 1824720846415195070;
table64[6] = 10146954863818469695;
table64[5] = 6540431468213822223;
table64[4] = 5852966179208367520;
table64[3] = 11377336531689676612;
table64[2] = 13423764633381537080;
table64[1] = 11121756678078876229;
table64[0] = 1;


table32 = (uint64_t *) malloc(sizeof(uint64_t)*32);
table32[31] = 670140591;
table32[30] = 2313048419;
table32[29] = 3291539972;
table32[28] = 2137603675;
table32[27] = 524737060;
table32[26] = 3629737640;
table32[25] = 1071065124;
table32[24] = 3907579083;
table32[23] = 2437735315;
table32[22] = 2257634306;
table32[21] = 2236890122;
table32[20] = 1940601056;
table32[19] = 3034121302;
table32[18] = 422655725;
table32[17] = 2737047479;
table32[16] = 483268876;
table32[15] = 598622898;
table32[14] = 4241522890;
table32[13] = 593004969;
table32[12] = 2255373408;
table32[11] = 4018602200;
table32[10] = 1762847826;
table32[9] = 3893040025;
table32[8] = 4001182008;
table32[7] = 3310684310;
table32[6] = 2730148548;
table32[5] = 4023817678;
table32[4] = 553251785;
table32[3] = 987466155;
table32[2] = 46175753;
table32[1] = 1631470458;
table32[0] = 1;

table16 = (uint64_t *) malloc(sizeof(uint64_t)*16);
table16[15] = 25555;
table16[14] = 22037;
table16[13] = 43066;
table16[12] = 21281;
table16[11] = 43747;
table16[10] = 13;
table16[9] = 92;
table16[8] = 4364;
table16[7] = 15473;
table16[6] = 47584;
table16[5] = 16018;
table16[4] = 61362;
table16[3] = 34554;
table16[2] = 37061;
table16[1] = 44235;
table16[0] = 1;

table8 = (uint64_t *) malloc(sizeof(uint64_t)*8);

table8[7] = 113;
table8[6] = 212;
table8[5] = 7;
table8[4] = 18;
table8[3] = 11;
table8[2] = 78;
table8[1] = 215;
table8[0] = 1;

table4 = (uint64_t *)malloc(sizeof(uint64_t)*4);

table4[3] = 15;
table4[2] = 5;
table4[1] = 7;
table4[0] = 1;

table2 = (uint64_t *)malloc(sizeof(uint64_t)*2);

table2[1] = 3;
table2[0] = 1;

}


uint64_t *getTable(int m) {
	if(m==64) {
		return table64;
	}

	else if(m==32) {
		return table32;
	}
	else if(m==16) {
		return table16;
	}
	else if(m==8) {
		return table8;
	}
	else if(m==4) {
		return table4;
	}
//this means m==2 (m= 1 will not produce a polynomial field)
	else {
		return table2;
	}
}

void freeTable() {
	free(table64);
	free(table32);
	free(table16);
	free(table8);
	free(table4);
	free(table2);

}

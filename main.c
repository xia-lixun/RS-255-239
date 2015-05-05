/* 
 * File:   main.c
 * Author: lixun
 *
 * Reed Solomon (255, 239) Encoder.
 * 
 * Created on August 16, 2014, 10:32 PM
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>
#include "ecc.h"



#define ML (239 + NPAR)
//N = 255
//K = 239
//t = (N - K) / 2 = 8
//m = 8



static uint8_t x[255];    //polynomial division working buffer [x^254 x^253 x^252 ... x^16   x^15 ... x^2 x^1 x^0]
static uint8_t i2d[255];  //map index form [alpha^0 .. alpha^254] to symbols (decimal)
static uint8_t d2i[256];  //map symbols to index form
//static uint8_t g[17] = {1, 59, 13, 104, 189, 68, 209, 30, 8, 163, 65, 41, 229, 98, 50, 36, 59};
static uint8_t g[17] = {1, 118, 52, 103, 31, 104, 126, 187, 232, 17, 56, 183, 49, 100, 81, 44, 79};
static uint8_t q;
static uint16_t ind, ind1, ind2;





void rs255239_init(void) {
    
    //generate index-to-decimal form lookup table
    // P(x) = x^8 + x^4 + x^3 + x^2 + 1
    uint16_t Px = 0x11D;
    uint16_t a  = 0x1;
    
    i2d[0] = (uint8_t)a;
    for(int i = 1; i < 255; i++) {
        a = a << 1;
        if(a & 0x100) {
            a = a ^ Px;
        }
        i2d[i] = (uint8_t)a;
    }
    
    
    //generate decimal-to-index form lookup table
    //note that i2d[i] belongs to 1..255
    //d2i[0] is a padding slot for operation reduction
    for(int i = 0; i < 255; i++) {
        d2i[i2d[i]] = (uint8_t)i;
    }
    
}




void rs255239_pdiv(void) {
    
    //polynomial division from kx^238 to kx^0
    for(int i = 0; i < 239; i++) {
        if(x[i] == 0x0) {
            continue;
        } else {
            ind1 = (uint16_t) d2i[x[i]];
            for(int j = 0; j < 17; j++) {
                //convert from decimal to index
                ind2 = (uint16_t) d2i[g[j]];    
                ind  = ind1 + ind2;
                ind  = ind % 255;
                q  = i2d[ind];
                x[i+j] ^= q;
            }
        }
    }    
}




//msg[] must be at least of length 239 bytes
//the polynomial remainder will be stored in x[239 .. 254] for [... x^15 ... x^2 x^1 x^0]
//append the remainder to msg[] bytes to form 255 byte encoded block
void rs255239_encode(uint8_t msg[]) {
    
    memcpy(x, msg, 239);
    memset(x+239, 0x0, 16);
    rs255239_pdiv();
}



//code[] must be at least of length 255 bytes
//syndrome will be stored in x[239 .. 254]
//correct code word results in all zeros
void rs255239_synd(uint8_t code[]) {
    
    memcpy(x, code, 255);
    rs255239_pdiv();
}




/*
 * 
 */
int main(int argc, char** argv) {

    rs255239_init();
    
//    for(int k = 0; k < 255; k++) {
//        printf("%4x  %d\n", i2d[k], i2d[k]);
//    }
//    printf("====\n");
//    for(int k = 1; k < 256; k++) {
//        printf("%4x  %d\n", d2i[k], d2i[k]);
//    }
    
    uint8_t msg[256];
    
    //encode
    for(int i = 0; i < 256; i++) {
        msg[i] = (uint8_t)i;
    }
    rs255239_encode(msg+1);
    //form the code word
    for(int k = 239; k < 255; k++) {
        msg[k+1] = x[k];
    }
    
    
    //print the parity
    for(int k = 239; k < 255; k++) {
        printf("%d ", x[k]);
    }
    printf("\n");    
    

    
    //calculate the syndrome, must be zeros
    rs255239_synd(msg+1);
    for(int k = 239; k < 255; k++) {
        printf("%d ", x[k]);
    }
    printf("\n");        
    
    
    
    
    
    
    
    //add error polynomial
    initialize_ecc ();
    int erasures[16];
    int nerasures = 0;    
    
    for(int loop = 0; loop < 16; loop++) {
    printf("======\n");
    for(int i = 0; i < 256; i++) msg[i] = (uint8_t)i;
        
        rs255239_encode(msg+1);
        //form the code word
        for(int k = 239; k < 255; k++) {
            msg[k+1] = x[k];
        }
        printf("++");     //parity
        for(int k=240; k<256; k++ ) printf("%d ", msg[k]);
        printf("++\n");     
        
        int errpos = (rand() % 255)+1;
        msg[errpos] ^= rand() % 256;
        printf("<");     //parity
        printf("%d -> %d", errpos, msg[errpos]);
        printf(">\n");             
        
        rs255239_synd(msg+1);  //syndrome
        printf("--");        
        for(int k = 239; k < 255; k++) printf("%d ", x[k]);
        printf("--\n");        
        
        
//    printf("--");        
//    encode_data(msg+1, 239, x);
//    for(int k=0; k<256; k++ ) printf("%d ", x[k]);
//    printf("--\n");        
        
        
    printf("\ndecode:");
    //test the decoding capability with 3rd party tools
    decode_data(msg+1, ML);
    print_syndrome();
    debug_check_syndrome();
    if (check_syndrome () != 0) {
        correct_errors_erasures (msg+1, ML, nerasures, erasures);
        printf("<");     //parity
        printf("%d -> %d", errpos, msg[errpos]);
        printf(">\n");
    }
    }
      
    return (EXIT_SUCCESS);
}



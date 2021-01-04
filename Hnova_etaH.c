#include "../binary_c.h"

double Hnova_etaH(struct star_t * accretor,
                  struct stardata_t * RESTRICT const stardata)
{
    //based on data from WANG2018 
    //note on variable naming: higher && lower in the 'fit' variable name refers to the mass of the WD curve 
    //(eg. a 1.38 Msun WD curve will be the 'higher', compared to a 1.35 Msun WD curve.)
    
    //returns: eta_H
    // Required units for inputs: M_sun and M_sun/yr (for M_WD and the accretion rate respectively)

    double Mdotlog10=log10(Mdot_net(accretor));
    double M_WD = accretor->mass;
    double fitlowereta=0, fithighereta=0, M_WDlower=0, M_WDhigher=0;
    double interpetaH;

//### INTERP MASS //##
    if(M_WD>=0.6 && M_WD<=0.7){
        
        //0.6 curves
        if(Mdotlog10>-8.097){
            fitlowereta= (1.26882419*pow(Mdotlog10,3) + 30.51329427*pow(Mdotlog10,2)
                          + 244.71761203*Mdotlog10 + 654.81045119);
        }
        else if(Mdotlog10>-10){
            fitlowereta=0.0658*Mdotlog10 + 0.7992;
        }
        else{
            fitlowereta= (0.0375*Mdotlog10 + 0.515);
        }

            
        //0.7 curves
        if(Mdotlog10>-8.222){
            fithighereta= (0.60632932*pow(Mdotlog10,3) + 14.74271288*pow(Mdotlog10,2) +
                          119.52930871*Mdotlog10 + 323.39248424);
        }
        else if(Mdotlog10>-10){
            fithighereta= 0.0719*Mdotlog10 + 0.845;
        }
        else{
            fithighereta= (0.0375*Mdotlog10 + 0.489);
        }
    
        M_WDlower=0.6;
        M_WDhigher=0.7;
    }
    else if(M_WD>=0.7 && M_WD<=0.8){
        //0.7 curves
        if(Mdotlog10>-8.222){
            fitlowereta= (0.60632932*pow(Mdotlog10,3) + 14.74271288*pow(Mdotlog10,2) +
                          119.52930871*Mdotlog10 + 323.39248424);
        }
        else if(Mdotlog10>-10){
            fitlowereta= 0.0719*Mdotlog10 + 0.845;
        }
        else{
            fitlowereta= (0.0375*Mdotlog10 + 0.489);
        }
            
            
        //0.8 curves
        if(Mdotlog10>-8.222){
            fithighereta=( -0.50065077*pow(Mdotlog10,4) - 14.76592690*pow(Mdotlog10,3)
                         - 162.29194136*pow(Mdotlog10,2) - 786.96965746*Mdotlog10 - 1418.36919563);
        }
        else if(Mdotlog10>-10){
            fithighereta=0.0444*Mdotlog10 + 0.5501;
        }
        else{
            fithighereta= 0.0375*Mdotlog10 + 0.472;
        }
        
        M_WDlower=0.7;
        M_WDhigher=0.8;
    }
    else if(M_WD>=0.8 && M_WD<=0.9){
        //0.8 curves
        if(Mdotlog10>-8.222){
            fitlowereta=( -0.50065077*pow(Mdotlog10,4) - 14.76592690*pow(Mdotlog10,3)
                         - 162.29194136*pow(Mdotlog10,2) - 786.96965746*Mdotlog10 - 1418.36919563);
        }
        else if(Mdotlog10>-10){
            fitlowereta=0.0444*Mdotlog10 + 0.5501;
        }
        else{
            fitlowereta= 0.0375*Mdotlog10 + 0.472;
        }
            
        //0.9 curves
        if(Mdotlog10>-8.222){
            fithighereta=(0.74958587*pow(Mdotlog10,4) + 23.11975024*pow(Mdotlog10,3) + 
                          267.36253222*pow(Mdotlog10,2) + 1374.05805223*Mdotlog10 + 2648.42148035);
        }
        else if(Mdotlog10>-10){
            fithighereta=0.0321*Mdotlog10 + 0.4261;
        }
        else{
            fithighereta=  0.0375*Mdotlog10 + 0.464;
        }
        
        M_WDlower=0.8;
        M_WDhigher=0.9;
    }
        
    else if(M_WD>=0.9 && M_WD<=1.0){
        //0.9 curves
        if(Mdotlog10>-8.222){
            fitlowereta=(0.74958587*pow(Mdotlog10,4) + 23.11975024*pow(Mdotlog10,3) + 
                          267.36253222*pow(Mdotlog10,2) + 1374.05805223*Mdotlog10 + 2648.42148035);
        }
        else if(Mdotlog10>-10){
            fitlowereta=0.0321*Mdotlog10 + 0.4261;
        }
        else{
            fitlowereta=  0.0375*Mdotlog10 + 0.464;
        }
            
        //1.0 curves
        if(Mdotlog10>-7.301){
            fithighereta=(18.23435276*pow(Mdotlog10,4) + 514.63820396*pow(Mdotlog10,3) +
                          5445.75095214*pow(Mdotlog10,2) + 25606.40507043*Mdotlog10 + 45143.59657352);
        }
        else if(Mdotlog10>-10){
            fithighereta= 0.061*Mdotlog10 + 0.6977;
        }

        else{
            fithighereta=0.0375*Mdotlog10 + 0.459;
        }


        
        M_WDlower=0.9;
        M_WDhigher=1.0;
    }
    
    else if(M_WD>=1.0 && M_WD<=1.1){
        //1.0 curves
        if(Mdotlog10>-7.301){
            fitlowereta=(18.23435276*pow(Mdotlog10,4) + 514.63820396*pow(Mdotlog10,3) +
                          5445.75095214*pow(Mdotlog10,2) + 25606.40507043*Mdotlog10 + 45143.59657352);
        }
        else if(Mdotlog10>-10){
            fitlowereta= 0.061*Mdotlog10 + 0.6977;
        }

        else{
            fitlowereta=0.0375*Mdotlog10 + 0.459;
        }
            
            
        //1.1 curves
        if(Mdotlog10>-6.481){
            fithighereta=29.464*Mdotlog10 + 191.59;
        }
        
        else if(Mdotlog10>-10){
            fithighereta=(0.01257896*pow(Mdotlog10,6) + 0.63924936*pow(Mdotlog10,5)
                          + 13.49899565*pow(Mdotlog10,4) + 151.61793668*pow(Mdotlog10,3) +
                          955.30681729*pow(Mdotlog10,2) + 3201.56125284*Mdotlog10 + 4458.98741717);
        }
        else{
            fithighereta=0.0375*Mdotlog10 + 0.478;
        }


        
        M_WDlower=1.0;
        M_WDhigher=1.1;
    }

    else if(M_WD>=1.1 && M_WD<=1.2){
        //1.1 curves
        if(Mdotlog10>-6.481){
            fitlowereta=29.464*Mdotlog10 + 191.59;
        }
        
        else if(Mdotlog10>-10){
            fitlowereta=(0.01257896*pow(Mdotlog10,6) + 0.63924936*pow(Mdotlog10,5)
                          + 13.49899565*pow(Mdotlog10,4) + 151.61793668*pow(Mdotlog10,3) +
                          955.30681729*pow(Mdotlog10,2) + 3201.56125284*Mdotlog10 + 4458.98741717);
        }
        else{
            fitlowereta=0.0375*Mdotlog10 + 0.478;
        }

            
        //1.2 curves
        if(Mdotlog10>-6.398){
            fithighereta= 38.404*Mdotlog10 + 246.09;
        }
        
        else if(Mdotlog10>-10){
            fithighereta=( 0.01863275*pow(Mdotlog10,4) + 0.63375932*pow(Mdotlog10,3) +
                          8.03922144*pow(Mdotlog10,2) + 45.11100840*Mdotlog10 + 94.67330838);
        }
        else{
            fithighereta= 0.0375*Mdotlog10 + 0.426;
        }

        M_WDlower=1.1;
        M_WDhigher=1.2;
    }
    else if(M_WD>=1.2 && M_WD<=1.3){
        //1.2 curves
        if(Mdotlog10>-6.398){
            fitlowereta= 38.404*Mdotlog10 + 246.09;
        }
        
        else if(Mdotlog10>-10){
            fitlowereta=(0.01863275*pow(Mdotlog10,4) + 0.63375932*pow(Mdotlog10,3) +
                          8.03922144*pow(Mdotlog10,2) + 45.11100840*Mdotlog10 + 94.67330838);
        }
        else{
            fitlowereta= 0.0375*Mdotlog10 + 0.426;
        }

            
        //1.3 curves
        if(Mdotlog10>-6.301){
            fithighereta= 63.255*Mdotlog10 + 399.02;
        }
        
        else if(Mdotlog10>-10){
            fithighereta=(0.00761654*pow(Mdotlog10,6) + 0.38074886*pow(Mdotlog10,5) +
                          7.90781049*pow(Mdotlog10,4) + 87.34785494*pow(Mdotlog10,3) +
                          541.22637679*pow(Mdotlog10,2) + 1783.79685277*Mdotlog10 + 2443.49600493);
        }
        else{
            fithighereta= 0.0375*Mdotlog10 + 0.442;
        }

        M_WDlower=1.2;
        M_WDhigher=1.3;
    }
    else if(M_WD>=1.3 && M_WD<=1.35){
        //1.3 curves
        if(Mdotlog10>-6.301){
            fitlowereta= 63.255*Mdotlog10 + 399.02;
        }
        
        else if(Mdotlog10>-10){
            fitlowereta=(0.00761654*pow(Mdotlog10,6) + 0.38074886*pow(Mdotlog10,5) +
                          7.90781049*pow(Mdotlog10,4) + 87.34785494*pow(Mdotlog10,3) +
                          541.22637679*pow(Mdotlog10,2) + 1783.79685277*Mdotlog10 + 2443.49600493);
        }
        else{
            fitlowereta= 0.0375*Mdotlog10 + 0.442;
        }

            
        //1.35 curves
        if(Mdotlog10>-6.301){
            fithighereta= 13.644*Mdotlog10 + 86.505;
        }
        
        else if(Mdotlog10>-10){
            fithighereta=( 0.02356744*pow(Mdotlog10,4) + 0.79605073*pow(Mdotlog10,3) +
                          10.02225648*pow(Mdotlog10,2) + 55.77497830*Mdotlog10 + 116.04245802);
        }
        else{
            fithighereta=  0.0375*Mdotlog10 + 0.511;
        }


        M_WDlower=1.3;
        M_WDhigher=1.35;
    }
//### EXTERAP MASS //##   

    else if(M_WD>1.35){
        //1.35 curves
        if(Mdotlog10>-6.301){
            fitlowereta= 13.644*Mdotlog10 + 86.505;
        }
        
        else if(Mdotlog10>-10){
            fitlowereta=( 0.02356744*pow(Mdotlog10,4) + 0.79605073*pow(Mdotlog10,3) +
                          10.02225648*pow(Mdotlog10,2) + 55.77497830*Mdotlog10 + 116.04245802);
        }
        else{
            fitlowereta=  0.0375*Mdotlog10 + 0.511;
        }
        
        
        //1.44 curves
        if(Mdotlog10>-6.301){
            fithighereta=13.644*Mdotlog10 + 86.545;
        }
        
        else if(Mdotlog10>-10){
            fithighereta=(0.02356744*pow(Mdotlog10,4) + 0.79605073*pow(Mdotlog10,3) +
                          10.02225648*pow(Mdotlog10,2) + 55.77497830*Mdotlog10 + 116.08245802);
        }
        
        else{
            fithighereta= 0.0375*Mdotlog10 + 0.551;
        }
        
        M_WDlower=1.35;
        M_WDhigher=1.44;
    }
    else if(M_WD<0.6){
       //0.2 curves
        if(Mdotlog10>-10.8){
            fitlowereta=(0.03375668*pow(Mdotlog10,6) + 1.96936922*pow(Mdotlog10,5) +
                         47.85435449*pow(Mdotlog10,4) + 619.95702757*pow(Mdotlog10,3) +
                         4516.25622905*pow(Mdotlog10,2) + 17541.05736755*Mdotlog10 + 28378.89454540);
        }
        else{
            fitlowereta=(0.0375*Mdotlog10 + 0.551);
        }
        
        
        //0.6 curves
        if(Mdotlog10>-8.097){
            fithighereta= (1.26882419*pow(Mdotlog10,3) + 30.51329427*pow(Mdotlog10,2)
                          + 244.71761203*Mdotlog10 + 654.81045119);
        }
        else if(Mdotlog10>-10){
            fithighereta=0.0658*Mdotlog10 + 0.7992;
        }
        else{
            fithighereta= (0.0375*Mdotlog10 + 0.515);
        }
        
        M_WDlower=0.2;
        M_WDhigher=0.6;
    }
//### RETURN && LIMIT //##
    // allow shifting of efficiency by amount nova_eta_shift
    // printf("before (fithighereta) H %g \n",
    //     fithighereta);

    fithighereta += stardata->preferences->nova_eta_shift;
    fitlowereta += stardata->preferences->nova_eta_shift;
    
    // printf("after (fithighereta) H %g \n",
    //     fithighereta);
    if(fithighereta>1){
        fithighereta=1;
        }
    
    if(fitlowereta>1){
        fitlowereta=1;
        }
    
    interpetaH=interp_lin(M_WD,M_WDlower,M_WDhigher,fitlowereta,fithighereta);

    return interpetaH;
}
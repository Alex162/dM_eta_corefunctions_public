#include "../binary_c.h"

double Henova_dMHe(struct star_t * accretor)
{
    //From formulae fitted to data from Kato et al 2018 (Figure 1,
    //data obtained from M. Kato, private communication) for M_WD>1.0.
    //From formulae fitted to data from Piersanti et al 2014
    //(supplementary data B7) for M_WD <1.0
    
    //note on variable naming: higher && lower in the 'fit' variable name refers to the mass of the WD curve 
    //(eg. a 1.38 Msun WD curve will be the 'higher', compared to a 1.35 Msun WD curve.)
    
    //returns: dMHe
    
    // Required units for inputs: M_sun and M_sun/yr (for M_WD and the accretion rate respectively)
    
    double Mdotlog10=log10(Mdot_gain(accretor));
    double M_WD = accretor->mass;
    
    double M_WD_Chand=1.44;///DO NOT CHANGE, this is just a parameter for the fits, NOT anything to do with binary_c
    double M_Heig_chand=log(8e-7);//He
    double fithigherMtranlog10=0, fitlowerMtranlog10=0, M_WDlower=0, M_WDhigher=0;
    
    //initialising dMHenova:
    double interpvalueMtran=0;
    
//#####Mdotlog10<-6 K2018: use quartic fits.###############
    if(M_WD<=1.38 && M_WD>=1.35){
        //1.38 curves
        if(Mdotlog10<-6){
            fithigherMtranlog10=(0.54437094*pow(Mdotlog10,4) + 14.68324636*pow(Mdotlog10,3) +
                                 148.24381872*pow(Mdotlog10,2) + 663.22092056*Mdotlog10 + 1103.57417218);
                                 }
        else{
            fithigherMtranlog10= -0.1056*Mdotlog10 - 5.6752;
        }
            
        //1.35 curves
        if(Mdotlog10<-6){
            fitlowerMtranlog10=(0.47596972*pow(Mdotlog10,4) + 12.61002433*pow(Mdotlog10,3)
                                + 125.08958464*pow(Mdotlog10,2) + 549.88858261*Mdotlog10 + 898.29715596);
        }
        else{
            fitlowerMtranlog10= -0.0934*Mdotlog10 - 5.3203;
        }
        
        M_WDhigher=1.38;
        M_WDlower=1.35;
    }


            
    else if(M_WD<=1.35 && M_WD>=1.3){
        //1.35 curves
        if(Mdotlog10<-6){
            fithigherMtranlog10=(0.47596972*pow(Mdotlog10,4) + 12.61002433*pow(Mdotlog10,3)
                                + 125.08958464*pow(Mdotlog10,2) + 549.88858261*Mdotlog10 + 898.29715596);
        }
        else{
            fithigherMtranlog10= -0.0934*Mdotlog10 - 5.3203;
        }
            
        //1.3 curves:
        if(Mdotlog10<-6){
            fitlowerMtranlog10=(1.16410162*pow(Mdotlog10,4) + 30.86542279*pow(Mdotlog10,3)
                                + 306.21947654*pow(Mdotlog10,2) + 1346.49353679*Mdotlog10 + 2208.95849065);
        }
        else{
            fitlowerMtranlog10= -0.1574*Mdotlog10 - 5.3083;
        }
        M_WDhigher=1.35;
        M_WDlower=1.3;
    }

        
    else if(M_WD<=1.30 && M_WD>=1.25){
        //1.3 curves
        if(Mdotlog10<-6){
            fithigherMtranlog10=(1.16410162*pow(Mdotlog10,4) + 30.86542279*pow(Mdotlog10,3)
                                + 306.21947654*pow(Mdotlog10,2) + 1346.49353679*Mdotlog10 + 2208.95849065);
        }
        else{
            fithigherMtranlog10= -0.1574*Mdotlog10 - 5.3083;
        }
            
        //1.25 curves
        if(Mdotlog10<-6){
            fitlowerMtranlog10=(1.68244304*pow(Mdotlog10,4) + 45.15739606*pow(Mdotlog10,3)
                                + 453.76725009*pow(Mdotlog10,2) + 2022.35971770*Mdotlog10 + 3368.09376231);
        }
        else{
            fitlowerMtranlog10= -0.1274*Mdotlog10 - 4.7732;
        }
        
        
        M_WDhigher=1.3;
        M_WDlower=1.25;
    }

            
            
    else if(M_WD<=1.25 && M_WD>=1.20){
        //1.25 curves
        if(Mdotlog10<-6){
            fithigherMtranlog10=(1.68244304*pow(Mdotlog10,4) + 45.15739606*pow(Mdotlog10,3)
                                + 453.76725009*pow(Mdotlog10,2) + 2022.35971770*Mdotlog10 + 3368.09376231);
        }
        else{
            fithigherMtranlog10= -0.1274*Mdotlog10 - 4.7732;
        };
        
        //1.20 curves
        if(Mdotlog10<-6){
            fitlowerMtranlog10= (0.91755524*pow(Mdotlog10,4) + 24.10100929*pow(Mdotlog10,3)
                                 + 236.98730163*pow(Mdotlog10,2) + 1033.14600808*Mdotlog10 + 1680.17435446);
        }
        else{
            fitlowerMtranlog10= -0.0542*Mdotlog10 - 4.1937;
        }
            
        M_WDhigher=1.25;
        M_WDlower=1.20;
    }
            
    else if(M_WD<=1.20 && M_WD>=1.10){
        //1.20 curves
        if(Mdotlog10<-6){
            fithigherMtranlog10= (0.91755524*pow(Mdotlog10,4) + 24.10100929*pow(Mdotlog10,3)
                                 + 236.98730163*pow(Mdotlog10,2) + 1033.14600808*Mdotlog10 + 1680.17435446);
        }
        else{
            fithigherMtranlog10= -0.0542*Mdotlog10 - 4.1937;
        }
        //1.10 curves
        if(Mdotlog10<-6){
            fitlowerMtranlog10= (1.72662751*pow(Mdotlog10,4) + 45.88762486*pow(Mdotlog10,3) 
                                 + 456.63345197*pow(Mdotlog10,2) + 2015.69281770*Mdotlog10 + 3325.89047580);
        }
        else{
            fitlowerMtranlog10= -0.0937*Mdotlog10 - 4.0218;
        }
        M_WDhigher=1.20;
        M_WDlower=1.10;
    }

    else if(M_WD<=1.10 && M_WD>=1.02){
//############1.10 FROM KATO, 1.02 FROM PIERSANTI. THIS IS THE BRIDGE############
        //1.10 curves
        if(Mdotlog10<-6){
            fithigherMtranlog10= (1.72662751*pow(Mdotlog10,4) + 45.88762486*pow(Mdotlog10,3) 
                                 + 456.63345197*pow(Mdotlog10,2) + 2015.69281770*Mdotlog10 + 3325.89047580);
        }
        else{
            fithigherMtranlog10= -0.0937*Mdotlog10 - 4.0218;
        }
        //1.02 curves
        if(Mdotlog10>-7.09){
            fitlowerMtranlog10= -0.9181*Mdotlog10 - 8.6514;
        }
        else{
            fitlowerMtranlog10=-2.0406*Mdotlog10 - 16.597;
        }
        
        M_WDhigher=1.10;
        M_WDlower=1.02;
    }
//############FROM HERE ON OUT ITS JUST PIERSANTI. All of these can be done adequately with linear fits################.
    else if(M_WD<=1.02 && M_WD>=0.92){
        //NOTE){ else if(MEANS THE ABOVE IS EFFECTIVELY <=-1, RATHER THAN <=1.02.
        //0.92 curves
        if(Mdotlog10>-7.30){
            fitlowerMtranlog10= -0.8228*Mdotlog10 - 7.8311;
        }
        else{
            fitlowerMtranlog10= -2.0714*Mdotlog10 - 16.929;
        }

        //1.02 curves
        if(Mdotlog10>-7.09){
            fithigherMtranlog10= -0.9181*Mdotlog10 - 8.6514;
        }
        else{
            fithigherMtranlog10=-2.0406*Mdotlog10 - 16.597;
        }
        
        M_WDlower=0.92;
        M_WDhigher=1.02;
    }
    
    else if(M_WD<=0.92 && M_WD>=0.81){
        //0.81 curves
        if(Mdotlog10>-7.40){
            fitlowerMtranlog10= -1.2462*Mdotlog10 - 10.729;
        }
        else{
            fitlowerMtranlog10=-2.7234*Mdotlog10 - 21.695;
        }
            
        //0.92 curves
        if(Mdotlog10>-7.30){
            fithigherMtranlog10= -0.8228*Mdotlog10 - 7.8311;
        }
        else{
            fithigherMtranlog10= -2.0714*Mdotlog10 - 16.929;
        }
        
        M_WDlower=0.81;
        M_WDhigher=0.92;
    }

            
    else if(M_WD<=0.81 && M_WD>=0.7){
        //0.7 curves
        if(Mdotlog10>-7.39){
            fitlowerMtranlog10 = -0.9928*Mdotlog10 - 8.8399;
        }
        else{
            fitlowerMtranlog10 = -2.8393*Mdotlog10 - 22.535;
        }

        //0.81 curves
        if(Mdotlog10>-7.40){
            fithigherMtranlog10= -1.2462*Mdotlog10 - 10.729;
        }
        else{
            fithigherMtranlog10=-2.7234*Mdotlog10 - 21.695;
        }
        
        M_WDlower=0.7;
        M_WDhigher=0.81;
    }
  
            
    else if(M_WD<=0.7 && M_WD>=0.6){
        //0.6curves
        if(Mdotlog10>-7.39){
            fitlowerMtranlog10= -1.2658*Mdotlog10 - 10.79;
        }
        else{
            fitlowerMtranlog10= -2.3314*Mdotlog10 - 18.668;
        }

        
        //0.7 curves
        if(Mdotlog10>-7.39){
            fithigherMtranlog10 = -0.9928*Mdotlog10 - 8.8399;
        }
        else{
            fithigherMtranlog10 = -2.8393*Mdotlog10 - 22.535;
        }
        
        M_WDlower=0.6;
        M_WDhigher=0.7;
    }


            
            
//####extrapolation in mass####){
    
    else if(M_WD<0.6){
        //EXTRAPOLATE based on 0.6 from P2014 && my phantom point at 0.2
        fitlowerMtranlog10=0.2;
        //0.6curves
        if(Mdotlog10>-7.39){
            fithigherMtranlog10= -1.2658*Mdotlog10 - 10.79;
        }
        else{
            fithigherMtranlog10= -2.3314*Mdotlog10 - 18.668;
        }
        
        M_WDlower=0.2;
        M_WDhigher=0.6;
    }

    //###Mass extrapolation past 1.38 with Mdot interpolation K2018){
    else if(M_WD>1.38){
        fithigherMtranlog10=M_Heig_chand;
        //1.38 curves
        if(Mdotlog10<-6){
            fitlowerMtranlog10=(0.54437094*pow(Mdotlog10,4) + 14.68324636*pow(Mdotlog10,3) +
                                 148.24381872*pow(Mdotlog10,2) + 663.22092056*Mdotlog10 + 1103.57417218);
        }
        else{
            fitlowerMtranlog10= -0.1056*Mdotlog10 - 5.6752;
        }
        
        M_WDhigher=M_WD_Chand;
        M_WDlower=1.38;
    }

        
    
//####interpolating && returning

    interpvalueMtran=pow(10,interp_lin(M_WD,M_WDlower,M_WDhigher,fitlowerMtranlog10,fithigherMtranlog10));

        
    return interpvalueMtran;
}
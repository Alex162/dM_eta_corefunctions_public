#include "../binary_c.h"

double Hedet_dMHe(struct star_t * accretor)
{
    //From formulae fitted to data from Kato et al 2018 (Figure 1,
    //data obtained from M. Kato, private communication) for M_WD>1.0.
    //From formulae fitted to data from Piersanti et al 2014
    //(supplementary data B7) for M_WD <1.0
    
    //note on variable naming: higher && lower in the 'fit' variable name refers to the mass of the WD curve 
    //(eg. a 1.38 Msun WD curve will be the 'higher', compared to a 1.35 Msun WD curve.)
    
    //returns: dMHe
    
    // Required units for inputs: M_sun and M_sun/yr (for M_WD and the accretion rate respectively)
    
    double Mdotlog10=log10(Mdot_net(accretor));
    double M_WD = accretor->mass;
    
    double M_WD_Chand=1.44;///DO NOT CHANGE, this is just a parameter for the fits, NOT anything to do with binary_c
    double M_Heig_chand=log(8e-7);//He
    double fithigherMtranlog10=0, fitlowerMtranlog10=0, M_WDlower=0, M_WDhigher=0;
    
    //initialising dMHedet:
    double interpvalueMtran=0;
////#Mdotlog10<-6 K2018:###############
    if(M_WD<=1.38 && M_WD>=1.35){
        //1.38 curves:
        if(Mdotlog10>-7.94){
            fithigherMtranlog10=(0.54437094*pow(Mdotlog10,4) + 14.68324636*pow(Mdotlog10,3) +
                                 148.24381872*pow(Mdotlog10,2) + 663.22092056*Mdotlog10 + 1103.57417218);
        }
        else{
            fithigherMtranlog10=-0.1433*Mdotlog10 - 4.079;
        }
        //1.35 curves:
        if(Mdotlog10>=-7.88){
            fitlowerMtranlog10=(0.47596972*pow(Mdotlog10,4) + 12.61002433*pow(Mdotlog10,3)
                                + 125.08958464*pow(Mdotlog10,2) + 549.88858261*Mdotlog10 + 898.29715596);
        }
        else{
            fitlowerMtranlog10= -0.0999*Mdotlog10 - 3.1982;
        }
        M_WDhigher=1.38;
        M_WDlower=1.35;
    }

            
    else if(M_WD<=1.35 && M_WD>=1.3){
        //1.35 curves:
        if(Mdotlog10>=-7.88){
            fithigherMtranlog10=(0.47596972*pow(Mdotlog10,4) + 12.61002433*pow(Mdotlog10,3)
                                + 125.08958464*pow(Mdotlog10,2) + 549.88858261*Mdotlog10 + 898.29715596);
        }
        else{
            fithigherMtranlog10= -0.0999*Mdotlog10 - 3.1982;
        }
            
        //1.30 curves:
        if(Mdotlog10>=-7.6920763){
            fitlowerMtranlog10=(1.16410162*pow(Mdotlog10,4) + 30.86542279*pow(Mdotlog10,3)
                                + 306.21947654*pow(Mdotlog10,2) + 1346.49353679*Mdotlog10 + 2208.95849065);
        }
        else{
            fitlowerMtranlog10=-0.0888*Mdotlog10 - 2.8977;
        }
        M_WDhigher=1.35;
        M_WDlower=1.3;
    }
        
    else if(M_WD<=1.30 && M_WD>=1.25){
        //1.30 curves:
        if(Mdotlog10>=-7.6920763){
            fithigherMtranlog10=(1.16410162*pow(Mdotlog10,4) + 30.86542279*pow(Mdotlog10,3)
                                + 306.21947654*pow(Mdotlog10,2) + 1346.49353679*Mdotlog10 + 2208.95849065);
        }
        else{
            fithigherMtranlog10=-0.0888*Mdotlog10 - 2.8977;
        }
        //1.25 curves
        if(Mdotlog10>=-7.62){
            fitlowerMtranlog10=(1.68244304*pow(Mdotlog10,4) + 45.15739606*pow(Mdotlog10,3)
                                + 453.76725009*pow(Mdotlog10,2) + 2022.35971770*Mdotlog10 + 3368.09376231);
        }
        else{
            fitlowerMtranlog10=-0.0957*Mdotlog10 - 2.8734;
        }

        M_WDhigher=1.3;
        M_WDlower=1.25;     
    }
            
            
    else if(M_WD<=1.25 && M_WD>=1.20){
        //1.25 curves
        if(Mdotlog10>=-7.62){
            fithigherMtranlog10=(1.68244304*pow(Mdotlog10,4) + 45.15739606*pow(Mdotlog10,3)
                                + 453.76725009*pow(Mdotlog10,2) + 2022.35971770*Mdotlog10 + 3368.09376231);
        }
        else{
            fithigherMtranlog10=-0.0957*Mdotlog10 - 2.8734;
        }
        
        //1.20 curves
        if(Mdotlog10>=-7.67){
            fitlowerMtranlog10= (0.91755524*pow(Mdotlog10,4) + 24.10100929*pow(Mdotlog10,3)
                                 + 236.98730163*pow(Mdotlog10,2) + 1033.14600808*Mdotlog10 + 1680.17435446);
        }
        else{
            fitlowerMtranlog10= -0.0646*Mdotlog10 - 2.1626;
        }

        M_WDhigher=1.25;
        M_WDlower=1.20;
    }

    else if(M_WD<=1.20 && M_WD>=1.10){
        //1.20 curves
        if(Mdotlog10>=-7.67){
            fithigherMtranlog10= (0.91755524*pow(Mdotlog10,4) + 24.10100929*pow(Mdotlog10,3)
                                 + 236.98730163*pow(Mdotlog10,2) + 1033.14600808*Mdotlog10 + 1680.17435446);
        }
        else{
            fithigherMtranlog10= -0.0646*Mdotlog10 - 2.1626;
        }
        //1.10 curves
        if(Mdotlog10>=-7.56){
            fitlowerMtranlog10= (1.72662751*pow(Mdotlog10,4) + 45.88762486*pow(Mdotlog10,3) 
                                 + 456.63345197*pow(Mdotlog10,2) + 2015.69281770*Mdotlog10 + 3325.89047580);
        }
        else{
            fitlowerMtranlog10=-0.0731*Mdotlog10 - 2.1163;
        }
        M_WDhigher=1.20;
        M_WDlower=1.10;
    }
        
//######1.10 FROM KATO, 1.02 FROM PIERSANTI. THIS IS THE BRIDGE.############
    else if(M_WD<=1.10 && M_WD>=1.02 && Mdotlog10<-6){

        //1.10 curves (K2018);
        if(Mdotlog10>=-7.56){
            fithigherMtranlog10= (1.72662751*pow(Mdotlog10,4) + 45.88762486*pow(Mdotlog10,3) 
                                 + 456.63345197*pow(Mdotlog10,2) + 2015.69281770*Mdotlog10 + 3325.89047580);
        }
        else{
            fithigherMtranlog10=-0.0731*Mdotlog10 - 2.1163;
        }
        
        //1.02 curves (P2014);
        if(Mdotlog10>=-7.30){
            fitlowerMtranlog10= -2.0406*Mdotlog10 - 16.597;
        }
        else if(Mdotlog10>=-8.82){
            fitlowerMtranlog10=(-0.52114281*pow(Mdotlog10,3) - 13.26176747*pow(Mdotlog10,2)
                                   - 112.71207541*Mdotlog10 - 320.54975268);
        }
        else{
            fitlowerMtranlog10=-0.0282*Mdotlog10 - 0.7734;
        }
            
            
        M_WDhigher=1.10;
        M_WDlower=1.02;
    }
        

        
//############PiersantiMdotlog10.
    else if(M_WD<=1.02 && M_WD>=0.92){
        
        //0.92 curves
        if(Mdotlog10>=-7.52){
            fitlowerMtranlog10= -2.0714*Mdotlog10 - 16.929;
        }
        else if(Mdotlog10>=-8.82){
            fitlowerMtranlog10= (-2.22457174*pow(Mdotlog10,4) - 74.06815881*pow(Mdotlog10,3) - 924.25967441*pow(Mdotlog10,2)
                                 - 5123.20522589*Mdotlog10 - 10644.68164309);
        }
        else{
            fitlowerMtranlog10= (-0.0238*Mdotlog10 - 0.6265);
        }        
                //1.02 curves (P2014);
        if(Mdotlog10>=-7.30){
            fithigherMtranlog10= -2.0406*Mdotlog10 - 16.597;
        }
        else if(Mdotlog10>=-8.82){
            fithigherMtranlog10=(-0.52114281*pow(Mdotlog10,3) - 13.26176747*pow(Mdotlog10,2)
                                   - 112.71207541*Mdotlog10 - 320.54975268);
        }
        else{
            fithigherMtranlog10=-0.0282*Mdotlog10 - 0.7734;
        }
            
        
        M_WDlower=0.92;
        M_WDhigher=1.02;
    }
    else if(M_WD<=0.92 && M_WD>=0.81){
        //0.81 curves
        if(Mdotlog10>=-7.60){
            fitlowerMtranlog10=-2.7234*Mdotlog10 - 21.695;
        }
        else if(Mdotlog10>=-8.82){
            fitlowerMtranlog10= (-0.89021703*pow(Mdotlog10,3) - 22.38723555*pow(Mdotlog10,2) 
                                 - 187.76315351*Mdotlog10 - 525.64755404);
        }
        else{
            fitlowerMtranlog10= (-0.0265*Mdotlog10 - 0.5781);
        }
        
        
        //0.92 curves
        if(Mdotlog10>=-7.52){
            fithigherMtranlog10= -2.0714*Mdotlog10 - 16.929;
        }
        else if(Mdotlog10>=-8.82){
            fithigherMtranlog10= (-2.22457174*pow(Mdotlog10,4) - 74.06815881*pow(Mdotlog10,3) - 924.25967441*pow(Mdotlog10,2)
                                 - 5123.20522589*Mdotlog10 - 10644.68164309);
        }
        else{
            fithigherMtranlog10= (-0.0238*Mdotlog10 - 0.6265);
        }
            
        M_WDlower=0.81;
        M_WDhigher=0.92;
    }
            
    else if(M_WD<=0.81 && M_WD>=0.7){
        //0.7 curves
        if(Mdotlog10>=-7.70){
            fitlowerMtranlog10=-2.8393*Mdotlog10 - 22.535;
        }
        else if(Mdotlog10>=-9){
            fitlowerMtranlog10 = (-0.68502223*pow(Mdotlog10,4) - 23.27090657*pow(Mdotlog10,3) - 296.28435839*pow(Mdotlog10,2)
                                  - 1675.83253716*Mdotlog10 - 3553.76035429);
        }
        else{
            fitlowerMtranlog10 = (-0.0251*Mdotlog10 - 0.4658);
        }
        
        
        //0.81 curves
        if(Mdotlog10>=-7.60){
            fithigherMtranlog10=-2.7234*Mdotlog10 - 21.695;
        }
        else if(Mdotlog10>=-8.82){
            fithigherMtranlog10= (-0.89021703*pow(Mdotlog10,3) - 22.38723555*pow(Mdotlog10,2) 
                                 - 187.76315351*Mdotlog10 - 525.64755404);
        }
        else{
            fithigherMtranlog10= (-0.0265*Mdotlog10 - 0.5781);
        }
            
        M_WDlower=0.7;
        M_WDhigher=0.81;
    }
            
    else if(M_WD<=0.7 && M_WD>=0.6){
        //0.6 curves
        if(Mdotlog10>=-7.82){
            fitlowerMtranlog10=-2.3314*Mdotlog10 - 18.668;
        }

        else if(Mdotlog10>=-9){
            fitlowerMtranlog10= (-0.13105692*pow(Mdotlog10,3) - 3.35859429*pow(Mdotlog10,2) 
                                 - 28.84546767*Mdotlog10 - 83.28360749);
        }
        else{
            fitlowerMtranlog10= -0.0303*Mdotlog10 - 0.4544;
        }
        
                //0.7 curves
        if(Mdotlog10>=-7.70){
            fithigherMtranlog10=-2.8393*Mdotlog10 - 22.535;
        }
        else if(Mdotlog10>=-9){
            fithigherMtranlog10 = (-0.68502223*pow(Mdotlog10,4) - 23.27090657*pow(Mdotlog10,3) - 296.28435839*pow(Mdotlog10,2)
                                  - 1675.83253716*Mdotlog10 - 3553.76035429);
        }
        else{
            fithigherMtranlog10 = (-0.0251*Mdotlog10 - 0.4658);
        }
            
        M_WDlower=0.6;
        M_WDhigher=0.7;
    }

            
            
    //####extrapolation in mass####){
    
    else if(M_WD<0.6){
        //EXTRAPOLATE based on 0.6 from P2014 && my phantom point at 0.2
        fitlowerMtranlog10=0.2;
    
        
        //0.6 curves
        if(Mdotlog10>=-7.82){
            fithigherMtranlog10=-2.3314*Mdotlog10 - 18.668;
        }

        else if(Mdotlog10>=-9){
            fithigherMtranlog10= (-0.13105692*pow(Mdotlog10,3) - 3.35859429*pow(Mdotlog10,2) 
                                 - 28.84546767*Mdotlog10 - 83.28360749);
        }
        else{
            fithigherMtranlog10= -0.0303*Mdotlog10 - 0.4544;
        }
        
        M_WDlower=0.2;
        M_WDhigher=0.6;
    }

    //###Mass extrapolation past 1.38 with Mdot interpolation K2018){
    else if(M_WD>1.38){
        fithigherMtranlog10=M_Heig_chand;
    
        
        //1.38 curves){
        if(Mdotlog10>-7.94){
            fitlowerMtranlog10=(0.54437094*pow(Mdotlog10,4) + 14.68324636*pow(Mdotlog10,3) +
                                 148.24381872*pow(Mdotlog10,2) + 663.22092056*Mdotlog10 + 1103.57417218);
        }
        else{
            fitlowerMtranlog10=-0.1433*Mdotlog10 - 4.079;
        }
        
        M_WDhigher=M_WD_Chand;
        M_WDlower=1.38;
    }

        
    
//####interpolating and returning
    
    interpvalueMtran=pow(10,interp_lin(M_WD,M_WDlower,M_WDhigher,fitlowerMtranlog10,fithigherMtranlog10));
        
    return interpvalueMtran;
}
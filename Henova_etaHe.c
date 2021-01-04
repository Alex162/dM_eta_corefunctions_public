#include "../binary_c.h"

double Henova_etaHe(struct star_t * accretor,
                    struct stardata_t * RESTRICT const stardata)
{
    //based on data from WU2017
    //note on variable naming: higher && lower in the 'fit' variable name refers to the mass of the WD curve 
    //(eg. a 1.38 Msun WD curve will be the 'higher', compared to a 1.35 Msun WD curve.)
    
    //returns: eta_He
    // Required units for inputs: M_sun and M_sun/yr (for M_WD and the accretion rate respectively)

    double Mdotlog10=log10(Mdot_net(accretor));
    double M_WD = accretor->mass;
    double fitlowereta=0, fithighereta=0, M_WDlower=0, M_WDhigher=0;
    double interpetaHe;

    
//### INTERP MASS //##

    
    if(M_WD>=0.7 && M_WD<=0.8){
        //0.7 curves
        if(Mdotlog10>-6.921){
            fitlowereta= (128.69653622*pow(Mdotlog10,4) + 3498.23532301*pow(Mdotlog10,3) +
                          35657.38069353*pow(Mdotlog10,2) + 161530.75025466*Mdotlog10 + 274398.30600877);
        }
        else{
            fitlowereta= 0.0271*Mdotlog10 + 0.4214;
        }
            
            
        //0.8 curves
        if(Mdotlog10>-6.495){
            fithighereta=(-2.81727160*pow(Mdotlog10,4) + 842.92305137*pow(Mdotlog10,3) +
                          17025.57779911*pow(Mdotlog10,2) + 111404.51221858*Mdotlog10 + 241315.71447017);
        }
        else{
            fithighereta= 0.0376*Mdotlog10 + 0.5774;
        }
        
        M_WDlower=0.7;
        M_WDhigher=0.8;
    }

    else if(M_WD>=0.8 && M_WD<=0.9){
        //0.8 curves
        if(Mdotlog10>-6.495){
            fitlowereta=(-2.81727160*pow(Mdotlog10,4) + 842.92305137*pow(Mdotlog10,3) +
                          17025.57779911*pow(Mdotlog10,2) + 111404.51221858*Mdotlog10 + 241315.71447017);
        }
        else{
            fitlowereta= 0.0376*Mdotlog10 + 0.5774;
        }
        
            
            
        //0.9 curves
        if(Mdotlog10>-6.585){
            fithighereta=(371.49410540*pow(Mdotlog10,4) + 9601.91634884*pow(Mdotlog10,3) +
                          93061.48887580*pow(Mdotlog10,2) + 400844.20977850*Mdotlog10 + 647425.88293061);
        }
        else{
            fithighereta=  0.0328*Mdotlog10 + 0.4598;
        }

        
        M_WDlower=0.8;
        M_WDhigher=0.9;
    }
        
    else if(M_WD>=0.9 && M_WD<=1.0){
        //0.9 curves
        if(Mdotlog10>-6.585){
            fitlowereta=(371.49410540*pow(Mdotlog10,4) + 9601.91634884*pow(Mdotlog10,3) +
                          93061.48887580*pow(Mdotlog10,2) + 400844.20977850*Mdotlog10 + 647425.88293061);
        }
        else{
            fitlowereta=  0.0328*Mdotlog10 + 0.4598;
        }
            
        //1.0 curves
        if(Mdotlog10>-6.155){
            fithighereta=(8.3115*Mdotlog10 + 51.674);
        }
        else if(Mdotlog10>-7.046){
            fithighereta= 0.42950525*pow(Mdotlog10,2) + 6.08106227*Mdotlog10 + 21.66204702;
        }
        else{
            fithighereta= 0.0274*Mdotlog10 + 0.331;
        }



        
        M_WDlower=0.9;
        M_WDhigher=1.0;
    }
    
    else if(M_WD>=1.0 && M_WD<=1.1){
        //1.0 curves
        if(Mdotlog10>-6.155){
            fitlowereta=(8.3115*Mdotlog10 + 51.674);
        }
        else if(Mdotlog10>-7.046){
            fitlowereta= 0.42950525*pow(Mdotlog10,2) + 6.08106227*Mdotlog10 + 21.66204702;
        }
        else{
            fitlowereta= 0.0274*Mdotlog10 + 0.331;
        }
            
            
        //1.1 curves
        if(Mdotlog10>-5.921){
            fithighereta=8.3856*Mdotlog10 + 50.358;
        }
        else if(Mdotlog10>-7.097){
            fithighereta=(-0.19075473*pow(Mdotlog10,4) - 4.92869779*pow(Mdotlog10,3) -
                          47.22636728*pow(Mdotlog10,2) - 198.28033421*Mdotlog10 - 306.30153647);
        }
        else{
            fithighereta=0.0275*Mdotlog10 + 0.2987;
        }


        
        M_WDlower=1.0;
        M_WDhigher=1.1;
    }
        
    else if(M_WD>=1.1 && M_WD<=1.2){
        //1.1 curves
        if(Mdotlog10>-5.921){
            fitlowereta=8.3856*Mdotlog10 + 50.358;
        }
        else if(Mdotlog10>-7.097){
            fitlowereta=(-0.19075473*pow(Mdotlog10,4) - 4.92869779*pow(Mdotlog10,3) -
                          47.22636728*pow(Mdotlog10,2) - 198.28033421*Mdotlog10 - 306.30153647);
        }
        else{
            fitlowereta=0.0275*Mdotlog10 + 0.2987;
        }

            
        //1.2 curves
        if(Mdotlog10>-5.854){
            fithighereta=  6.7216*Mdotlog10 + 40.146;
        }
        else if(Mdotlog10>-7.222){
            fithighereta=( 0.61101079*pow(Mdotlog10,4) + 16.16260323*pow(Mdotlog10,3) +
                          160.39740390*pow(Mdotlog10,2) + 708.15701427*Mdotlog10 + 1174.46712731);
        }
        else{
            fithighereta=  0.0278*Mdotlog10 + 0.2994;
        }

        M_WDlower=1.1;
        M_WDhigher=1.2;
    }

    else if(M_WD>=1.2 && M_WD<=1.25){
        //1.2 curves
        if(Mdotlog10>-5.854){
            fitlowereta=  6.7216*Mdotlog10 + 40.146;
        }
        else if(Mdotlog10>-7.222){
            fitlowereta=( 0.61101079*pow(Mdotlog10,4) + 16.16260323*pow(Mdotlog10,3) +
                          160.39740390*pow(Mdotlog10,2) + 708.15701427*Mdotlog10 + 1174.46712731);
        }
        else{
            fitlowereta=  0.0278*Mdotlog10 + 0.2994;
        }

            
        //1.25 curves
        if(Mdotlog10>-5.824){
            fithighereta= 7.3403*Mdotlog10 + 43.604;
        }
        
        else if(Mdotlog10>-7.222){
            fithighereta=(0.70445443*pow(Mdotlog10,4) + 18.61070131*pow(Mdotlog10,3) +
                          184.34484390*pow(Mdotlog10,2) + 811.82891151*Mdotlog10 + 1342.09332478);
        }
        else{
            fithighereta=  0.0278*Mdotlog10 + 0.2924;
        }

        M_WDlower=1.2;
        M_WDhigher=1.25;
    }
    
    else if(M_WD>=1.25 && M_WD<=1.3){
        //1.25 curves
        if(Mdotlog10>-5.824){
            fitlowereta= 7.3403*Mdotlog10 + 43.604;
        }
        
        else if(Mdotlog10>-7.222){
            fitlowereta=(0.70445443*pow(Mdotlog10,4) + 18.61070131*pow(Mdotlog10,3) +
                          184.34484390*pow(Mdotlog10,2) + 811.82891151*Mdotlog10 + 1342.09332478);
        }
        else{
            fitlowereta=  0.0278*Mdotlog10 + 0.2924;
        }
            
            
        //1.3 curves
        if(Mdotlog10>-5.824){
            fithighereta=  5.2446*Mdotlog10 + 31.397;
        }
        else if(Mdotlog10>-7.222){
            fithighereta=( 0.35511452*pow(Mdotlog10,4) + 9.43073062*pow(Mdotlog10,3) +
                          93.97568200*pow(Mdotlog10,2) + 416.89622203*Mdotlog10 + 695.71187519);
        }
        else{
            fithighereta=  0.0278*Mdotlog10 + 0.2754;
        }


        M_WDlower=1.25;
        M_WDhigher=1.3;
    }
        
    else if(M_WD>=1.3 && M_WD<=1.35){
        //1.3 curves
        if(Mdotlog10>-5.824){
            fitlowereta=  5.2446*Mdotlog10 + 31.397;
        }
        else if(Mdotlog10>-7.222){
            fitlowereta=( 0.35511452*pow(Mdotlog10,4) + 9.43073062*pow(Mdotlog10,3) +
                          93.97568200*pow(Mdotlog10,2) + 416.89622203*Mdotlog10 + 695.71187519);
        }
        else{
            fitlowereta=  0.0278*Mdotlog10 + 0.2754;
        }
            
            
        //1.35 curves
        if(Mdotlog10>-5.824){
            fithighereta= 4.742*Mdotlog10 + 28.497;
        }

        else if(Mdotlog10>-7.301){
            fithighereta=(0.17534165*pow(Mdotlog10,2) + 2.81395280*Mdotlog10 + 11.31042978);
        }
        else{
            fithighereta=0.028*Mdotlog10 + 0.315;
        }



        M_WDlower=1.3;
        M_WDhigher=1.35;
    }

//### EXTERAP MASS //##   

    else if(M_WD>1.35){
        //1.35 curves
        if(Mdotlog10>-5.824){
            fitlowereta= 4.742*Mdotlog10 + 28.497;
        }
        else if(Mdotlog10>-7.301){
            fitlowereta=(0.17534165*pow(Mdotlog10,2) + 2.81395280*Mdotlog10 + 11.31042978);
        }
        else{
            fitlowereta=0.028*Mdotlog10 + 0.315;
        }
        
        
        //1.44 curves
        if(Mdotlog10>-5.824){
            fithighereta=4.742*Mdotlog10 + 28.597;
        }
        
        else if(Mdotlog10>-7.301){
            fithighereta=(0.17534165*pow(Mdotlog10,2) + 2.81395280*Mdotlog10 + 11.41042978);
        }
        
        else{
            fithighereta=  0.028*Mdotlog10 + 0.415;
        }

        
        M_WDlower=1.35;
        M_WDhigher=1.44;
    }
        
    else if(M_WD<0.7){
       //0.2 curves
        if(Mdotlog10>-7.421){
            fitlowereta=(128.69653622*pow(Mdotlog10,4) + 3755.62839546*pow(Mdotlog10,3) +
                         41097.77848236*pow(Mdotlog10,2) + 199876.15570843*Mdotlog10 + 364523.39925814);
        }
        else{
            fitlowereta=( 0.0284*Mdotlog10 + 0.4944);
        }
        
        
        //0.7 curves
        if(Mdotlog10>-6.921){
            fithighereta= (128.69653622*pow(Mdotlog10,4) + 3498.23532301*pow(Mdotlog10,3) +
                          35657.38069353*pow(Mdotlog10,2) + 161530.75025466*Mdotlog10 + 274398.30600877);
        }
        else{
            fithighereta= 0.0271*Mdotlog10 + 0.4214;
        }
        
        M_WDlower=0.2;
        M_WDhigher=0.7;
    }

//### RETURN && LIMIT //##
    // allow shifting of efficiency by amount nova_eta_shift
    fithighereta += stardata->preferences->nova_eta_shift;
    fitlowereta += stardata->preferences->nova_eta_shift;
    if(fithighereta>1){
        fithighereta=1;
        }
    
    if(fitlowereta>1){
        fitlowereta=1;
        }
    
    interpetaHe=interp_lin(M_WD,M_WDlower,M_WDhigher,fitlowereta,fithighereta);

    return interpetaHe;
}
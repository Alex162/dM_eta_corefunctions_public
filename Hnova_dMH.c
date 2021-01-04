#include "../binary_c.h"

double Hnova_dMH(struct star_t * accretor)
{
    /*
    from formulae fitted to data from Kato et al. 2014
    #Returns: dMH
    
    Required units for inputs: M_sun and M_sun/yr (for M_WD and the accretion rate respectively.)
    */
    double log10Mdot=log10(Mdot_gain(accretor));
    double M_WD = accretor->mass;
    //my ASSUMPTION for the limiting Mchand=1.44 case:
    //ignition mass is 5e-8 Msun and independent of Mdot.
    double M_WD_Chand=1.44;
    double M_Hig_chand=5e-8; //
    double dMH=0;
    double M_WD_eneg4=0, M_WD_6eneg5=0, M_WD_3eneg5=0, M_WD_1p7eneg5=0;
    double M_WD_1eneg5=0, M_WD_6eneg6=0, M_WD_3eneg6=0, M_WD_1eneg6=0, M_WD_3eneg7=0;
    if(log10Mdot>=-7.7){
        M_WD_eneg4= -0.01136821*pow(log10Mdot,2) - 0.28779422*log10Mdot - 0.87275137;
        M_WD_6eneg5=-0.00223874*pow(log10Mdot,3) - 0.07980979*pow(log10Mdot,2) - 0.95727391*log10Mdot - 2.88226600;
        M_WD_3eneg5=-0.02513162*pow(log10Mdot,2) - 0.51596425*log10Mdot - 1.55913869;
        M_WD_1p7eneg5= -0.00286937*pow(log10Mdot,3) - 0.09804557*pow(log10Mdot,2) - 1.11466878*log10Mdot - 3.04208220;
        M_WD_1eneg5=-0.00271579*pow(log10Mdot,3) - 0.09167551*pow(log10Mdot,2) - 1.02778587*log10Mdot - 2.59280133;
        M_WD_6eneg6= -0.00281513*pow(log10Mdot,3) - 0.09252102*pow(log10Mdot,2) - 1.00663332*log10Mdot - 2.35355189;
        M_WD_3eneg6=-0.00162294*pow(log10Mdot,3) - 0.06120773*pow(log10Mdot,2) - 0.71424226*log10Mdot - 1.34053269;
        M_WD_1eneg6= -0.00284309*pow(log10Mdot,3) - 0.08134364*pow(log10Mdot,2) - 0.77925207*log10Mdot - 1.13067755;
        M_WD_3eneg7= -0.01995254*pow(log10Mdot,2) - 0.31057708*log10Mdot + 0.17515477;
    }
    else if(log10Mdot>=-9){
        M_WD_eneg4=0.03123628*pow(log10Mdot,3) + 0.75350110*pow(log10Mdot,2) + 5.94211003*log10Mdot + 16.00580843;
        M_WD_6eneg5=-0.00223874*pow(log10Mdot,3) - 0.07980979*pow(log10Mdot,2) - 0.95727391*log10Mdot - 2.88226600;
        M_WD_3eneg5=0.03007799*pow(log10Mdot,3) + 0.69444233*pow(log10Mdot,2) + 5.20316206*log10Mdot + 13.54164455;
        M_WD_1p7eneg5= -0.00286937*pow(log10Mdot,3) - 0.09804557*pow(log10Mdot,2) - 1.11466878*log10Mdot - 3.04208220;
        M_WD_1eneg5=-0.00271579*pow(log10Mdot,3) - 0.09167551*pow(log10Mdot,2) - 1.02778587*log10Mdot - 2.59280133;
        M_WD_6eneg6= -0.00281513*pow(log10Mdot,3) - 0.09252102*pow(log10Mdot,2) - 1.00663332*log10Mdot - 2.35355189;
        M_WD_3eneg6=-0.00162294*pow(log10Mdot,3) - 0.06120773*pow(log10Mdot,2) - 0.71424226*log10Mdot - 1.34053269;
        M_WD_1eneg6= -0.00284309*pow(log10Mdot,3) - 0.08134364*pow(log10Mdot,2) - 0.77925207*log10Mdot - 1.13067755;
        M_WD_3eneg7=-0.0062 * log10Mdot + 1.3336;
    }
    else if(log10Mdot<-9){
        M_WD_eneg4=-0.0225 * log10Mdot + 0.5856;
        M_WD_6eneg5=-0.0246 * log10Mdot + 0.6755;
        M_WD_3eneg5=-0.0224 * log10Mdot + 0.8371;
        M_WD_1p7eneg5=-0.0173 * log10Mdot + 0.9779;
        M_WD_1eneg5=-0.0148 * log10Mdot + 1.0727;
        M_WD_6eneg6=-0.0141 * log10Mdot + 1.1325;
        M_WD_3eneg6=-0.0092 * log10Mdot + 1.2307;
        M_WD_1eneg6= -0.0047 * log10Mdot + 1.3239;
        M_WD_3eneg7=-0.0062 * log10Mdot + 1.3336;
    }
    
    //interpolation
    if((M_WD<=M_WD_6eneg5) && (M_WD>=M_WD_eneg4)){
        // printf("dMH is between 6e-5 and e-4 Msun ");
        dMH=interp_lin(M_WD,M_WD_eneg4,M_WD_6eneg5,1e-4,6e-5);
    }
    else if((M_WD<=M_WD_3eneg5) && (M_WD>=M_WD_6eneg5)){
        // printf("M_WD is between 3e-5 and 6e-5 Msun ");
        dMH=interp_lin(M_WD,M_WD_6eneg5,M_WD_3eneg5,6e-5,3e-5);
    }    
    else if((M_WD<=M_WD_1p7eneg5) && (M_WD>=M_WD_3eneg5)){
        // printf("M_WD is between 1.7e-5 and 3e-5 Msun ");
        dMH=interp_lin(M_WD,M_WD_3eneg5,M_WD_1p7eneg5,3e-5,1.7e-5);
    }    
    else if((M_WD<=M_WD_1eneg5) && (M_WD>=M_WD_1p7eneg5)){
        // printf("M_WD is between 1e-5 and 1.7e-5 Msun ");
        dMH=interp_lin(M_WD,M_WD_1p7eneg5,M_WD_1eneg5,1.7e-5,1e-5);
    }    
    else if((M_WD<=M_WD_6eneg6) && (M_WD>=M_WD_1eneg5)){
        // printf("M_WD is between 6e-6 and 1e-5 Msun ");
        dMH=interp_lin(M_WD,M_WD_1eneg5,M_WD_6eneg6,1e-5,6e-6);
    }    
    else if((M_WD<=M_WD_3eneg6) && (M_WD>=M_WD_6eneg6)){
        // printf("M_WD is between 3e-6 and 6e-6 Msun ");
        dMH=interp_lin(M_WD,M_WD_6eneg6,M_WD_3eneg6,6e-6,3e-6);
    }
    else if((M_WD<=M_WD_1eneg6) && (M_WD>=M_WD_3eneg6)){
        // printf("M_WD is between 1e-6 and 3e-6 Msun ");
        dMH=interp_lin(M_WD,M_WD_3eneg6,M_WD_1eneg6,3e-6,1e-6);
    }
    else if((M_WD<=M_WD_3eneg7) && (M_WD>=M_WD_1eneg6)){
        // printf("M_WD is between 3e-7 and 1e-6 Msun ");
        dMH=interp_lin(M_WD,M_WD_1eneg6,M_WD_3eneg7,1e-6,3e-7);
    }
    
        // Some EXPRAPOLATION cases
    else if(M_WD<M_WD_eneg4){
        //dMH is greater than e-4 Msun... use my assumption/approximation
        //for limiting behaviour approaching M_WD=0.2.
        //dMH=interp_lin(M_WD,M_WD_eneg4,M_WD_6eneg5,1e-4,6e-5)
        dMH=interp_lin(M_WD,0.2,M_WD_eneg4,4e-4,1e-4);
    }
    else if(M_WD>M_WD_3eneg7){
        //dMH is less than 3e-7 Msun... use my
        //assumption/approximation for limiting behaviour on Mchand
        dMH=interp_lin(M_WD,M_WD_3eneg7,M_WD_Chand,3e-7,M_Hig_chand);
    }
    // printf("%g\n",dMH);
    return dMH;
}
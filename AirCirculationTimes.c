#include "udf.h"
#include "sg_mphase.h"
#include "cxndsearch.h"
#include <stdio.h>

#define NORMAL 1.6
#define CHANGE 3.2
#define FLUID_INDEX 9
#define H2_INDEX 0
#define CH4_INDEX 1
#define DECTION_X 100
#define DECTION_Y 2.9
#define LIMIT 0.0004

int flag = 0;

DEFINE_PROFILE(AirCirculation,thread,index)
{
	face_t f;
	Domain*domain= Get_Domain(1);  
	Material *mix_mat = mixture_material(domain); 
	Thread*tc=Lookup_Thread(domain, FLUID_INDEX);
	cell_t c;
    cell_t jiance;
	CX_Cell_Id *cx_cell;
	real molofch4 = 0;
	real ch4_v=0;
	real total_v=0;
	real ch4_vof = 0;
	
	double Pt_to_find[2]={DECTION_X,DECTION_Y};
	ND_Search*domain_table = NULL;
	domain_table = CX_Start_ND_Point_Search(domain_table,TRUE,-1); 
	cx_cell = CX_Find_Cell_With_Point(domain_table, Pt_to_find, 0);
    jiance = RP_CELL(cx_cell);
	Message0("jiance=%d\n",jiance);
	
	
	begin_c_loop(c,tc)
	{
		if(jiance==c){
		int i;
		real Mw_CH4;
        real c_v;
		real total_mole=0;
		Material *spe_mat = NULL;
		real CH4_mole_fract = 0;
		mixture_species_loop(mix_mat, spe_mat,i) 
		{
			real Mwi=MATERIAL_PROP(spe_mat,PROP_mwi); 
			if(i==CH4_INDEX) Mw_CH4=Mwi;  
			total_mole += C_YI(c,tc,i)/Mwi; 
		}
		CH4_mole_fract=(C_YI(c,tc,CH4_INDEX)/Mw_CH4)/total_mole; 
		c_v = C_VOLUME(c,tc);
		ch4_v += CH4_mole_fract*c_v;
		total_v += c_v;
		}
	}
	end_c_loop(c,tc);
	ch4_vof  = ch4_v/total_v;

	Message("mole of ch4 on jiance=%g\n",ch4_vof);

	if(ch4_vof>LIMIT){
		flag=1;
	}
	Message("flag=%d\n",flag);
    begin_f_loop(f,thread)
    {
        if(flag==0)
        {
            F_PROFILE(f,thread,index)=NORMAL;
        }
        else
        {
            F_PROFILE(f,thread,index)=CHANGE;
        }
    }
    end_f_loop(f,thread)
}
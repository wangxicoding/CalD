#include "Element.h"
#include "math.h"
using namespace std;

double slqh(double x, double y)//矢量求和函数
{
	double xy;
	xy=sqrt(x*x+y*y);
	return xy;
}

double sl_angel(double x, double y)//根据向量求角度，x，y分别表示XY方向的分量
{
	double angel;
	if(x==0)
	{
		if(y>0)
		{
			angel=0.5*PI;//Y轴正方向
		}
		else
		{
			angel=1.5*PI;//Y轴负方向
		}
	}    //过滤分母为零的情况，以免出现不必要的错误
	else
	{
		if(x>0)
		{
			if(y>0)
			{
				angel=fabs(atan(y/x));//第一象限角
			}
			else
			{
				angel=2*PI-fabs(atan(fabs(y)/x));//第四象限角
			}
		}
		else
		{
			if(y>0)
			{
				angel=PI-fabs(atan(y/fabs(x)));//第二象限角
			}
			else
			{
				angel=PI+fabs(atan(fabs(y)/fabs(x)));//第三象限角
			}

		}
		
	}
	if(fabs(angel-2*PI)<0.00001)//过滤360度情况
		{
			angel=0;
		}
	return angel;//返回向量角度，单位为弧度制
}


double CElement::AnsysElement()
{
    m_area=0.5*PI*m_r*m_r;
    m_moment=0.5*m_mass*m_r*m_r;
}
double CElement::Initial()
{
	//单元性质
	m_number=0; //单元编号初始化
	m_x=m_y=m_r=0.0; //单元x方向的坐标，单元y方向的坐标，单元的半径初始化
	m_velx=m_vely=m_rotate=0.0; //单元x，y方向的速度及单元转动速度初始化
	vt_t2[4]=vtt2[4]=mtt2=mt_t2=0.0; //vt_t2和mt_t2分别表示t-△t/2时刻的速度和转动速度，vtt2和mtt2分别表示t+△t/2时刻的速度和转动速度初始化
	m_mass=m_moment=0.0; //单元的质量和转动惯量初始化
	m_disx=m_disy=m_diso=0.0; //分别表示单元x，y方向的位移增量及转角增量初始化
	m_disn=m_diss=0.0; //分别表示单元法向方向位移增量，切向方向位移增量初始化
	m_ag=0.0; //地震加速度记录，从地震数据中读取初始化
	m_Feqx=m_Feqy=0.0; //x方向和y方向地震力初始化
	m_Fxsum=m_Fysum=m_Msum=0.0; //x和y方向和合力及合力矩初始化
	//材料性质参数
	m_young=0.0;
	m_possion=0.0;
	m_thickness=0.0;
	m_weight=0.0;
}

void CElement::calContactForce(CElement* p1, CElement* p2)//计算接触力有弹簧，p2单元作用到p1单元上的力，同时也储存到p2的接触力容器中，便于计算p2的接触力
{
   
	
	int i; //循环变量
	double l_xy; //代表两颗粒之间绝对距离
	double l_xc; //l_xc表示接触点到p1单元圆心的距离
	double lx, ly, c_x, c_y; //
	double sl_u[3]; //定义接触力矢量数组
	double sumR; //两个单元半径和
	double cos_theta,sin_theta; //定义两个单元的圆心连线与x轴正方向的余弦值和正弦值
	double delta_n,delta_s; //法向位移增量，切向位移增量
	double b1, b2;//法向力增量，切向力增量，法向刚度增量，切向刚度增量

	sl_u[0]=(p2->m_x)-(p1->m_x);
	sl_u[1]=(p2->m_y)-(p1->m_y);

	l_xy=slqh(sl_u[0],sl_u[1]);//调用求和程序，给定两个单元的x，y坐标值，计算两个单元的绝对距离
	cos_theta=(p2->m_x-p1->m_x)/l_xy;//计算两单元型心连线与x轴的夹角余弦值
	sin_theta=(p2->m_y-p1->m_y)/l_xy;//计算两单元型心连线与x轴的夹角正弦值

        l_xc=0.5*(p1->m_r-p2->m_r+l_xy); //

	lx=l_xc*cos_theta; //
	ly=l_xc*sin_theta; //

	c_x=lx+p1->m_x; //接触点c的x坐标
	c_y=ly+p1->m_y; //接触点c的y坐标

	b1=p2->m_disx-p1->m_disx; //
	b2=p2->m_disy-p1->m_disy; //
	//对应第2页第4步
	p1->m_disn=b1*cos_theta+b2*sin_theta；//计算法向方向位移增量
	p1->m_diss=-b1*sin_theta+b2*cos_theta-((p1->m_r)*(p1->m_diso)+(p2->m_r)*(p2->m_diso))；//切向方向位移增量
	//对应第2页第5步
	p1->p_contract.m_fn=p1->p_contract.m_fn-kn*p1->m_disn;//求法向力fn
	p1->p_contract.m_fs=p1->p_contract.m_fn+ks*p1->m_diss;//求切向力fs
	//对应第2页第6步
	p1->p_contract.m_dn=-beta*kn*(p1->m_disn);//求法向阻尼力Dn，beta为参数，设置为0.6
	p1->p_contract.m_ds=-beta*ks*(p1->m_diss);//求切向阻尼力Ds
	//对应第3页第7步
	p1->tau=C+mu*(p1->p_contract.m_fn);//求极限剪力值tau，C为常数,这个tau用于后面的判断，在DiscreteElement.cpp中
	//求合力及合力矩
	p1->m_Fxsum=-(p1->p_contract.m_fs+p1->p_contract.m_ds)*sin_theta-(p1->p_contract.m_fn+p1->p_contract.m_dn)*cos_theta;
        p1->m_Fysum=(p1->p_contract.m_fs+p1->p_contract.m_ds)*cos_theta-(p1->p_contract.m_fn+p1->p_contract.m_dn)*sin_theta;
	p1->m_Msum=p1->m_Fysum*(c_x-p1->m_x)-p1->m_Fxsum*(c_y-p1->m_y);
	//作用到p2单元上的力为
	p2->m_Fxsum=-p1->m_Fxsum;
	p2->m_Fysum=-p1->m_Fysum;
	p2->m_Msum=-p1->m_Msum;
}
void CElement::calCollisionForce(CElement* p1, CElement* p2) //结算碰撞力，无弹簧
{
	int i; //循环变量
	double l_xy; //代表两颗粒之间绝对距离
	double sl_u[3]; //定义接触力矢量数组
	double cos_theta,sin_theta; //定义两个单元的圆心连线与x轴正方向的余弦值和正弦值
	double delta_n,delta_s; //法向位移增量，切向位移增量
	double sumR; //两个单元半径和
	double e[2]; t[2];
	double f_n, f_s, d_n, d_s;//法向力增量，切向力增量，法向刚度增量，切向刚度增量
	double delta_Fx, delta_Fy;

	sumR=p1->m_r+p2->m_r; //求两个单元半径和
	
	sl_u[0]=(p2->m_x)-(p1->m_x);
	sl_u[1]=(p2->m_y)-(p1->m_y);

	l_xy=slqh(sl_u[0],sl_u[1]);//调用求和程序，给定两个单元的x，y坐标值，计算两个单元的绝对距离
	cos_theta=(p2->m_x-p1->m_x)/l_xy;//计算两单元型心连线与x轴的夹角余弦值
	sin_theta=(p2->m_y-p1->m_y)/l_xy;//计算两单元型心连线与x轴的夹角正弦值
	
	
	if(l_xy>sumR)
	{
		continue;
	}
	else
	{
		e[0]=((p2->m_x)-(p1->m_x))/l_xy; //cos_theta
		e[1]=((p2->m_y)-(p1->m_y))/l_xy; //sin_theta
		t[0]=e[1]; //sin_theta
		t[1]=-e[0]; //-cos_theta
	    
        for(i=0;i<2;i++)
	    {
		    vn_XD[0]=((p1->m_velx)-(p2->m_velx))*e[0];//单元相对速度x方向分量
            vn_XD[1]=((p1->m_vely)-(p2->m_vely))*e[1];//单元相对速度y方向分量
		    vs_XD[0]=(((p1->m_velx)-(p2->m_velx))*t[0])-(p1->m_rotate*p1->m_r+p2->m_rotate*p2->m_r));//相对速度
		    vs_XD[1]=(((p1->m_vely)-(p2->m_vely))*t[1])-(p1->m_rotate*p1->m_r+p2->m_rotate*p2->m_r));//相对速度
		    drt_n[i]=vn_XD[i]*delta;//相对位移增量
		    drt_s[i]=vs_XD[i]*delta;//相对位移增量
	    }
            vn_XD[3]=slqh(vn_XD[0],vn_XD[1]);//法向相对速度大小
            vs_XD[3]=slqh(vs_XD[0],vs_XD[1]);//切向相对速度大小
	        drt_n[3]=slqh(drt_n[0],drt_n[1]);//法向的位移增量
	        drt_s[3]=slqh(drt_s[0],drt_s[1]);//切向位移增量值

		    f_n=kn*drt_n[3];//法向力增量
		    f_s=ks*drt_s[3];//切向力增量
		    d_n=-beta*kn*vn_XD[3];//法向刚度增量
		    d_s=-beta*kn*vs_XD[3];//切向刚度增量

		    p1->p_contract.m_fn=p1->p_contract.m_fn+f_n+d_n;
		    p1->p_contract.m_fs=p1->p_contract.m_fs+f_s+d_s;
		    p1->tau=(p1->p_contract.m_fn)*mu+C;//库伦摩擦力

		    if(abs(p1->p_contract.m_fs)>abs(p1->tau))
		    {
			    p1->p_contract.m_fs=abs(p1->tau)*(p1->p_contract->m_fs/abs(p1->p_contract->m_fs));//方向与原摩阻力方向一致
		    }
		    else
			{
		        delta_Fx=p1->p_contract.m_fn*e[0]+p1->p_contract.m_fs*e[1];
		        delta_Fy=p1->p_contract.m_fn*e[1]-p1->p_contract.m_fs*e[0];
		        p1->p_contract.m_Fx=p1->p_contract->m_Fx+delta_Fx;//求x方向合力，未考虑外力和重力
				p1->p_contract.m_Fy=p1->p_contract->m_Fy+delta_Fy;//求x方向合力，未考虑外力和重力
		        p1->p_contract.m_M=p1->p_contract.m_M+p1->p_contract.m_fs*p1->m_r;//求合力矩
			}
	}
}

void CElement::cal_vtt2(CElement* p)
{
	double c1, c2;
	c1=1-alpha*drt/2;
	c2=1+alpha*drt/2;
	p->vtt2[0]=(p->vt_t2[0]*c1+p->m_Fxsum*drt/p->m_mass)/c2;
	p->vtt2[1]=(p->vt_t2[1]*c1+(p->m_Fysum/p->m_mass+g)*drt)/c2;
	p->mtt2=(p->mt_t2*c1+(p->m_Msum/p->moment)*drt)/c2;
        p->vtt2[2]=sl_angel(p->vtt2[0],p->vtt2[1]);
	p->vtt2[3]=slqh(p->vtt2[0],p->vtt2[1]);
}

void CElement::cal_vt(CElement* p)
{
	p->m_velx=(p->vtt2[0]+p->vt_t2[0])/2;
	p->m_vely=(p->vtt2[1]+p->vt_t2[1])/2;
	p->m_rotate=(p->mtt2+p->mt_t2)/2;
}

void CElement::cal_dis(CElement* p)
{
	p->m_disx=p->vtt2[0]*drt;
	p->m_disy=p->vtt2[1]*drt;
	p->m_diso=p->mtt2*drt;
}

void CElement::cal_utt(CElement* p)
{
	p->m_x=p->m_x+p->m_disx;
	p->m_y=p->m_y+p->m_disy;
	p->m_o=p->m_o+p->m_diso;
}

void CElement::change_of_data(CElement* p)
{
	int i;
	for(i=0; i<4; i++)
	{
		p->vt_t2[i]=p->vtt2[i];
	}
}

void CElement::union_lisan(CElement* p)
{
	cal_vtt2(p);
	cal_vt(p);
	cal_dis(p);
	cal_utt(p);
	change_of_data(p);
}

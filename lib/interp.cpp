#include "interp.h"

namespace rockseis {
// constructor
template<typename T>
Interp<T>::Interp()
{
    mode = LINEAR;
    pad = LN_PAD;
}

template<typename T>
Interp<T>::Interp(rs_interpmode _mode)
{
    mode = _mode;
    switch(mode){
        case LINEAR:
            pad = LN_PAD;
            break;
        case BSPLINE:
            pad = BS_PAD;
            break;
        case SINC:
            pad = SN_PAD;
            break;
    }
}

// destructor
template<typename T>
Interp<T>::~Interp(){
/* Do nothing*/
}

template<typename T>
T Interp<T>::interp_ln_int1d(T *mod, T d1)
/*<Linear interpolation function in 1d >*/
{
    T interpoint;

    interpoint=((1.0-d1)*mod[0] + d1*mod[1]);
    return (interpoint); 
}

template<typename T>
T Interp<T>::interp_bs_int1d(T *mod, T d1)
/*<B-spline interpolation function in 1d >*/
{
    T x_plus_1, one_minus_x;
    T a0, a1, a2, a3;
    T interpoint;

    x_plus_1=d1+1.0;
    one_minus_x=1.0-d1; 
    a0= ((-(1.0/6.0)*x_plus_1 +1.0)*x_plus_1 -2.0) * x_plus_1 + (4.0/3.0);
    a1=(0.5*d1 -1.0)*d1*d1 + (2.0/3.0);
    a2=(0.5*one_minus_x -1.0)*one_minus_x*one_minus_x + (2.0/3.0);
    a3=1.0 -a0-a1-a2;

    interpoint=a0*mod[0]+a1*mod[1]+a2*mod[2]+a3*mod[3];
    return (interpoint);
}

template<typename T>
T Interp<T>::sinc(const T x)
/*< Sinc function >*/
{
    if(x==0.0){
        return (1.0);
    }
    return (sin(PI*x)/(PI*x));
}

template<typename T>
T Interp<T>::interp_sinc_int1d(T *mod, T d1)
/*<Sinc interpolation function in 1d >*/
{
    int i;
    T interpoint;

    interpoint=0.0;
    for(i=0; i < 6 ; i++){
        interpoint=interpoint+mod[i]*sinc(((2.0-i)+d1))*sinc(((2.0-i)+d1)/3.0);
    }
    return (interpoint);
}


template<typename T>
void Interp<T>::interp(std::shared_ptr<Data2D<T>> from, std::shared_ptr<Data2D<T>> to)
{
    size_t n1, n2, ntr;
    T d1, d2;
    T o1, o2;
    T *data1 = from->getData();
    T *data2 = to->getData();
    Point2D<T> *Scoords1 = (from->getGeom())->getScoords();
    Point2D<T> *Gcoords1 = (from->getGeom())->getGcoords();
    Point2D<T> *Scoords2 = (to->getGeom())->getScoords();
    Point2D<T> *Gcoords2 = (to->getGeom())->getGcoords();
    n1 = from->getNt();
    n2 = to->getNt();
    ntr = from->getNtrace();
    if(ntr != to->getNtrace()) rs_error("Interp::interp: Number of traces in the two datasets mismatch."); 

    Index ind1(n1,ntr);
    Index ind2(n2,ntr);

    d1 = from->getDt();
    d2 = to->getDt();

    o1 = from->getOt();
    o2 = to->getOt();
    T val;
    long int i0;
    T x0;
    for (size_t i = 0; i < ntr; i++)
    {
        //Copy coordinates
        Scoords2[i].x = Scoords1[i].x;
        Scoords2[i].y = Scoords1[i].y;
        Gcoords2[i].x = Gcoords1[i].x;
        Gcoords2[i].y = Gcoords1[i].y;
        for(size_t j = 0; j < n2; j++)
        {
            val = o2 + j*d2;
            i0 = (long int) ((val-o1)/d1) - (this->pad - 1);
            if(i0 < 0){
                data2[ind2(j,i)] = 0.0;
                continue;
            }
            if(i0 > (n1 + 2*this->pad - 6)){
                data2[ind2(j,i)] = 0.0;
                continue;
            }
            x0 = (o1 + (i0+(this->pad - 1))*d1);
            switch(this->mode){
                case LINEAR:
                    data2[ind2(j,i)] = interp_ln_int1d(&data1[ind1(i0, i)], (val - x0)/d1);
                    break;
                case BSPLINE:
                    data2[ind2(j,i)] = interp_bs_int1d(&data1[ind1(i0, i)], (val - x0)/d1);
                    break;
                case SINC:
                    data2[ind2(j,i)] = interp_sinc_int1d(&data1[ind1(i0, i)], (val - x0)/d1);
                    break;
            }
        }
    }
}

template<typename T>
void Interp<T>::interp(std::shared_ptr<Data3D<T>> from, std::shared_ptr<Data3D<T>> to)
{
    size_t n1, n2, ntr;
    T d1, d2;
    T o1, o2;
    T *data1 = from->getData();
    T *data2 = to->getData();
    Point3D<T> *Scoords1 = (from->getGeom())->getScoords();
    Point3D<T> *Gcoords1 = (from->getGeom())->getGcoords();
    Point3D<T> *Scoords2 = (to->getGeom())->getScoords();
    Point3D<T> *Gcoords2 = (to->getGeom())->getGcoords();
    n1 = from->getNt();
    n2 = to->getNt();
    ntr = from->getNtrace();
    if(ntr != to->getNtrace()) rs_error("Interp::interp: Number of traces in the two datasets mismatch."); 

    Index ind1(n1,ntr);
    Index ind2(n2,ntr);

    d1 = from->getDt();
    d2 = to->getDt();

    o1 = from->getOt();
    o2 = to->getOt();
    T val;
    long int i0;
    T x0;
    for (size_t i = 0; i < ntr; i++)
    {
        //Copy coordinates
        Scoords2[i].x = Scoords1[i].x;
        Scoords2[i].y = Scoords1[i].y;
        Scoords2[i].z = Scoords1[i].z;
        Gcoords2[i].x = Gcoords1[i].x;
        Gcoords2[i].y = Gcoords1[i].y;
        Gcoords2[i].z = Gcoords1[i].z;
        for(size_t j = 0; j < n2; j++)
        {
            val = o2 + j*d2;
            i0 = (long int) ((val-o1)/d1) - (this->pad - 1);
            if(i0 < 0){
                data2[ind2(j,i)] = 0.0;
                continue;
            }
            if(i0 > (n1 + 2*this->pad - 1)){
                data2[ind2(j,i)] = 0.0;
                continue;
            }
            x0 = (o1 + (i0+(this->pad - 1))*d1);
            switch(this->mode){
                case LINEAR:
                    data2[ind2(j,i)] = interp_ln_int1d(&data1[ind1(i0, i)], (val - x0)/d1);
                    break;
                case BSPLINE:
                    data2[ind2(j,i)] = interp_bs_int1d(&data1[ind1(i0, i)], (val - x0)/d1);
                    break;
                case SINC:
                    data2[ind2(j,i)] = interp_sinc_int1d(&data1[ind1(i0, i)], (val - x0)/d1);
                    break;
            }
        }
    }
}


// =============== INITIALIZING TEMPLATE STRUCTS =============== //
template class Interp<float>;
template class Interp<double>;
}



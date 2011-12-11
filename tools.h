#ifndef TOOLS_H
#define TOOLS_H 
#include<iomanip>
#include<string>
#include<tr1/random>
using namespace std;
/*
template <class T>
string stringify(T x, int width=15, const char ch=' ')
 {
   std::ostringstream o;
   if (!(o << setw(width)<<setfill(ch)<<x))
     cerr<<"Bad coversion to string"<<endl;
   return o.str();
 }

*/

extern std::tr1::normal_distribution<double> normal;
extern std::tr1::ranlux64_base_01 eng;

template<class T>
T rand_vec(){
	T v;
	for(size_t i=0; i<T::dim; i++){
	 v(i)=normal(eng);
	}
	v.normalize();
	return v;
}

template<class T>
void mysort(const vector<T> &eval, size_t ind[] , size_t n){ 
 
        bool swapped; 
 
	for(size_t i=1; i<n; i++){
	ind[i]=i;
	}
	ind[0]=0;
        double val1, val2; 
        do{ 
                swapped = false; 
                for(size_t i=1; i<n-1; ++i){ 
                         val1=eval[ind[i]]->print_width;
                         val2=eval[ind[i+1]]->print_width;
                      if (val1> val2){ 
                        swap( ind[i], ind[i+1] ); 
                        swapped = true; 
                        } 
                        } 
                --n; 
        }while (swapped); 

}

class COffset
	{
	public:
	COffset(double _t=0, double _b=0):top(_t),bottom(_b){}
	double width(){return top+bottom;}
	double top, bottom;
 	private:
	};

template<class T>
class CLink
	{
	public:
	CLink(T *_h, int d=1):head(_h),length(d){ }
	int length;
	T *head;
 	private:
	};

#endif /* TOOLS_H */

#ifndef TOOLS_H
#define TOOLS_H 
#include<iomanip>
#include<string>
using namespace std;

template <class T>
string stringify(T x, int width=15, const char ch=' ')
 {
   std::ostringstream o;
   if (!(o << setw(width)<<setfill(ch)<<x))
     cerr<<"Bad coversion to string"<<endl;
   return o.str();
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

#endif /* TOOLS_H */

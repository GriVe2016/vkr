#ifndef SeparatorsH
#define SeparatorsH

#include "property.h"
//#################################################################

class TSeparator
{
    public:
     TSeparator(const long double lStartValue, const long double lEndValue);
     virtual ~TSeparator() {}
     virtual TSeparator &Copy() const = 0;

	 long double GetStartValue() const {return fStartValue;}
     void SetStartValue(const long double& N) {
         throw "Not implemented";
     }
     Property<TSeparator,
         long double,
         &TSeparator::GetStartValue,
         &TSeparator::SetStartValue> StartValue;

	 long double GetEndValue() const {return fEndValue;}
     void SetEndValue(const long double& N) {
         throw "Not implemented";
     }
     Property<TSeparator,
         long double,
         &TSeparator::GetEndValue,
         &TSeparator::SetEndValue> EndValue;

	 int GetEndIndex() const {return fEndIndex;}
     void SetEndIndex(const int& EndIndex) {
         this->fEndIndex = EndIndex;
     }
     Property<TSeparator,
         int,
         &TSeparator::GetEndIndex,
         &TSeparator::SetEndIndex> EndIndex;

	 int GetCountNodes() const {return fEndIndex+1;}
     void SetCountNodex(const int& CountNodes) {
         throw "Not implemented";
     }
     Property<TSeparator,
         int,
         &TSeparator::GetCountNodes,
         &TSeparator::SetCountNodex> CountNodes;
     
     const long double *GetDimension() const {return fDimension;}
     void SetDimension(const long double * const & Dimension) {
         throw "Not implemented";
     }
     Property<TSeparator,
         const long double *,
         &TSeparator::GetDimension,
         &TSeparator::SetDimension> Dimension;

     const long double *GetSeparation() const {return fSeparation;}
     void SetSeparation(const long double * const & Separation) {
         throw "Not implemented";
     }
     Property<TSeparator,
         const long double *,
         &TSeparator::GetSeparation,
         &TSeparator::SetSeparation> Separation;

    protected:
       const long double fStartValue; //Начало отрезка. Совпадает с первым узлом.
       const long double fEndValue; //Конец отрезка. Совпадает с последним узлом.
       int fEndIndex; //Номер последнего узла. Номер первого узла - 0.
       long double *fDimension; //Узлы разбиения.
       long double *fSeparation; //Части разбиения.

    private:
      TSeparator(const TSeparator &):fStartValue(0),fEndValue(0),
           StartValue(this), EndValue(this), EndIndex(this),
           CountNodes(this), Dimension(this), Separation(this) {}
      TSeparator &operator=(const TSeparator &) {return *this;}      
};

//#################################################################

class TUniformSeparator: public TSeparator
{
    public:
     TUniformSeparator(const long double lStartValue, const long double lEndValue, const long double lSeparationValue);
     virtual ~TUniformSeparator();
     virtual TUniformSeparator &Copy() const;

     long double GetSeparationValue() const {return fSeparationValue;}
     void SetSeparationValue(const long double& SeparationValue) {
         this->fSeparationValue = SeparationValue;
     }
     Property<TUniformSeparator,
         long double,
         &TUniformSeparator::GetSeparationValue,
         &TUniformSeparator::SetSeparationValue> SeparationValue;

    protected:
     long double fSeparationValue; //Шаг разбиения.
    private:
      TUniformSeparator(const TUniformSeparator &): TSeparator(0,0), fSeparationValue(0), SeparationValue(this) {}
      TUniformSeparator &operator=(const TUniformSeparator &) {return *this;}
};

//#################################################################

class TRandomSeparator: public TSeparator
{
    public:
     TRandomSeparator(int lCountNodes, const long double *lDimension);
     virtual ~TRandomSeparator();
     virtual TRandomSeparator &Copy() const;
    private:
      TRandomSeparator(const TRandomSeparator &): TSeparator(0,0) {}
      TRandomSeparator &operator=(const TRandomSeparator &) {return *this;}
};

//#################################################################

#endif

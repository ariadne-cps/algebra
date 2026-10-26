/***************************************************************************
 *            chebyshev_polynomial.tpl.hpp
 *
 *  Copyright  2008-20  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "numeric/numeric.hpp"

#include <iomanip>
#include <limits>

#include "numeric/rounding.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/expansion.hpp"
#include "algebra/series.hpp"
#include "algebra/differential.hpp"
#include "algebra/chebyshev_polynomial.hpp"
#include "utility/exceptions.hpp"

#include "algebra/expansion.inl.hpp"
#include "algebra/expansion.tpl.hpp"
#include "algebra/evaluate.tpl.hpp"
#include "algebra/algebra_operations.tpl.hpp"

#include "algebra/multi_index-noaliasing.hpp"
#include "algebra/vector.hpp"

namespace Ariadne {

template<class X> auto UnivariateChebyshevPolynomial<X>::constant(X const& x) -> ChebyshevPolynomial<X> {
    UnivariateChebyshevPolynomial<X> cmr(x.precision());
    cmr._terms.append(0u,x);
    return cmr;
}

template<class X> auto UnivariateChebyshevPolynomial<X>::constant(Number<P> c, PR pr) -> ChebyshevPolynomial<X> {
    ChebyshevPolynomial<X> cmr(pr);
    cmr._terms.append(0u,X(c,pr));
    return cmr;
}

template<class X> auto UnivariateChebyshevPolynomial<X>::coordinate(PR pr) -> ChebyshevPolynomial<X> {
    ChebyshevPolynomial<X> cmr(pr);
    cmr._terms.append(1u,X(1,pr));
    return cmr;
}

template<class X> auto UnivariateChebyshevPolynomial<X>::basis(SizeType k, PR pr) -> ChebyshevPolynomial<X> {
    ChebyshevPolynomial<X> cmr(pr);
    cmr._terms.append(static_cast<DegreeType>(k),X(1,pr));
    return cmr;
}


template<class X> auto UnivariateChebyshevPolynomial<X>::create_constant(Number<P> c) -> ChebyshevPolynomial<X> {
    return UnivariateChebyshevPolynomial<X>::constant(c,this->precision());
}


template<class X> auto
UnivariateChebyshevPolynomial<X>::sup_norm() const -> MagType<X> {
    MagType<X> r = mag(this->zero_coefficient());
    for(auto term : this->_terms) { r+=mag(term.coefficient()); }
    return r;
}


template<class X> auto
UnivariateChebyshevPolynomial<X>::apply(Pos, ChebyshevPolynomial<X> cm) -> ChebyshevPolynomial<X> {
    for(auto term : cm._terms) { term.coefficient() = +term.coefficient(); }
    return cm;
}

template<class X> auto
UnivariateChebyshevPolynomial<X>::apply(Neg, ChebyshevPolynomial<X> cm) -> ChebyshevPolynomial<X> {
    for(auto term : cm._terms) { term.coefficient() = -term.coefficient(); }
    return cm;
}

template<class X> auto
UnivariateChebyshevPolynomial<X>::apply(Add, ChebyshevPolynomial<X> cm, X const& s) -> ChebyshevPolynomial<X> {
    if (cm._terms.empty()) { cm._terms.append(0u,s); }
    else if (cm._terms.front().index()==0) { cm._terms.front().coefficient()+=s; }
    else { cm._terms.prepend(0u,s); }
    return cm;
}

template<class X> auto
UnivariateChebyshevPolynomial<X>::apply(Mul, ChebyshevPolynomial<X> cm, X const& s) -> ChebyshevPolynomial<X> {
    for(auto term : cm._terms) { term.coefficient() *= s; }
    return cm;
}

template<class X> auto
UnivariateChebyshevPolynomial<X>::apply(Sqr, ChebyshevPolynomial<X> cm) -> ChebyshevPolynomial<X> {
    return cm*cm;
}

/*
template<class X> template<class Y> requires Convertible<SumType<X,Y>,X> auto
UnivariateChebyshevPolynomial<X>::_add(ChebyshevPolynomial<X> cm, Y const& s) -> ChebyshevPolynomial<X> {
    if (cm._terms.empty()) { cm._terms.append(0u,cm._terms.zero_coefficient()+s); }
    else if (cm._terms.front().index()==0) { cm._terms.front().coefficient()+=s; }
    else { assert(false); }
    return cm;
}

template<class X> template<class Y> requires Convertible<DifferenceType<X,Y>,X> auto
UnivariateChebyshevPolynomial<X>::_sub(ChebyshevPolynomial<X> cm, Y const& s) -> ChebyshevPolynomial<X> {
    if (cm._terms.empty()) { cm._terms.append(0u,cm._terms.zero_coefficient()-s); }
    else if (cm._terms.front().index()==0) { cm._terms.front().coefficient()-=s; }
    else { assert(false); }
    return cm;
}

template<class X> template<class Y> requires Convertible<ProductType<X,Y>,X> auto
UnivariateChebyshevPolynomial<X>::_mul(ChebyshevPolynomial<X> cm, Y const& s) -> ChebyshevPolynomial<X> {
    for(auto term : cm._terms) { term.coefficient() *= s; }
    return cm;
}

template<class X> template<class Y> requires Convertible<QuotientType<X,Y>,X> auto
UnivariateChebyshevPolynomial<X>::_div(ChebyshevPolynomial<X> cm, Y const& s) -> ChebyshevPolynomial<X> {
    for(auto term : cm._terms) { term.coefficient() = term.coefficient()/s; }
    return cm;
}
*/

template<class X> auto
UnivariateChebyshevPolynomial<X>::apply(Add, ChebyshevPolynomial<X> const& cm1, ChebyshevPolynomial<X> const& cm2) -> ChebyshevPolynomial<X> {
    UnivariateChebyshevPolynomial<X> cmr(min(cm1.precision(),cm2.precision()));
    auto iter1=cm1._terms.begin();
    auto iter2=cm2._terms.begin();
    while(iter1!=cm1._terms.end() && iter2!=cm2._terms.end()) {
        if(iter1->index()==iter2->index()) {
            cmr._terms.append(iter1->index(),iter1->coefficient()+iter2->coefficient());
            ++iter1; ++iter2;
        } else if(iter1->index()<iter2->index()) {
            cmr._terms.append(iter1->index(),iter1->coefficient());
            ++iter1;
        } else { // iter1->index()>iter2->index()
            cmr._terms.append(iter2->index(),iter2->coefficient());
            ++iter2;
        }
    }
    while(iter1!=cm1._terms.end()) {
        cmr._terms.append(iter1->index(),iter1->coefficient());
        ++iter1;
    }
    while(iter2!=cm2._terms.end()) {
        cmr._terms.append(iter2->index(),iter2->coefficient());
        ++iter2;
    }
    return cmr;
}

template<class X> auto
UnivariateChebyshevPolynomial<X>::apply(Sub, ChebyshevPolynomial<X> const& cm1, ChebyshevPolynomial<X> const& cm2) -> ChebyshevPolynomial<X> {
    return add(cm1,neg(cm2));
}

template<class X> auto
UnivariateChebyshevPolynomial<X>::apply(Mul, ChebyshevPolynomial<X> const& cm1, ChebyshevPolynomial<X> const& cm2) -> ChebyshevPolynomial<X> {
    PR pr=min(cm1.precision(),cm2.precision());
    X z(pr);
    Array<X> buf(cm1._terms.back().index()+cm2._terms.back().index()+1u,z);
    for (auto term1 : cm1._terms) {
        for (auto term2 : cm2._terms) {
            X coefficient_product=term1.coefficient()*term2.coefficient();
            DegreeType i1=term1.index(); DegreeType i2=term2.index();
            buf[static_cast<SizeType>(i1>=i2 ? i1-i2 : i2-i1)]+=coefficient_product;
            buf[term1.index()+term2.index()]+=coefficient_product;
        }
    }
    UnivariateChebyshevPolynomial<X> cmr(pr);
    for (SizeType i=0; i!=buf.size(); ++i) {
        cmr._terms.append(static_cast<DegreeType>(i),hlf(buf[i]));
    }
    return cmr;
}

template<class X> auto
UnivariateChebyshevPolynomial<X>::operator() (X const& x) const -> X {
    return _evaluate(*this,x);
}

template<class X> auto
UnivariateChebyshevPolynomial<X>::_write(OutputStream& os) const -> OutputStream& {
    for (auto term : this->_terms) {
        String c = to_str(term.coefficient());
        if (c.find('-')>=c.find_first_of("0123456789")) { os << "+"; }
        os << c << "*T" << term.index();
    }
    return os;
}





template<class X> auto MultivariateChebyshevPolynomial<X>::constant(SizeType as, X const& c) -> ChebyshevPolynomial<X> {
    ChebyshevPolynomial<X> cmr(as,c.precision());
    cmr._terms.append(MultiIndex::zero(as),c);
    return cmr;
}

template<class X> auto MultivariateChebyshevPolynomial<X>::constant(SizeType as, Number<P> c, PR pr) -> ChebyshevPolynomial<X> {
    ChebyshevPolynomial<X> cmr(as,pr);
    cmr._terms.append(MultiIndex::zero(as),X(c,pr));
    return cmr;
}

template<class X> auto MultivariateChebyshevPolynomial<X>::coordinate(SizeType as, SizeType i, PR pr) -> ChebyshevPolynomial<X> {
    ChebyshevPolynomial<X> cmr(as,pr);
    cmr._terms.append(MultiIndex::unit(as,i),X(1,pr));
    return cmr;
}

template<class X> auto MultivariateChebyshevPolynomial<X>::basis(SizeType as, SizeType i, SizeType k, PR pr) -> ChebyshevPolynomial<X> {
    ChebyshevPolynomial<X> cmr(as,pr);
    cmr._terms.append(MultiIndex::unit(as,i)*static_cast<DegreeType>(k),X(1,pr));
    return cmr;
}


template<class X> auto MultivariateChebyshevPolynomial<X>::create_constant(Number<P> c) -> ChebyshevPolynomial<X> {
    return ChebyshevPolynomial<X>::constant(this->argument_size(),c,this->precision());
}


template<class X> auto
MultivariateChebyshevPolynomial<X>::sup_norm() const -> MagType<X> {
    MagType<X> r = mag(this->zero_coefficient());
    for(auto term : this->_terms) { r+=mag(term.coefficient()); }
    return r;
}

template<class X> auto
MultivariateChebyshevPolynomial<X>::apply(Pos, ChebyshevPolynomial<X> cm) -> ChebyshevPolynomial<X> {
    for(auto term : cm._terms) { term.coefficient() = +term.coefficient(); }
    return cm;
}

template<class X> auto
MultivariateChebyshevPolynomial<X>::apply(Neg, ChebyshevPolynomial<X> cm) -> ChebyshevPolynomial<X> {
    for(auto term : cm._terms) { term.coefficient() = -term.coefficient(); }
    return cm;
}

template<class X> auto
MultivariateChebyshevPolynomial<X>::apply(Add, ChebyshevPolynomial<X> cm, Scalar<X> const& s) -> ChebyshevPolynomial<X> {
    if (cm._terms.empty()) { cm._terms.append(MultiIndex::zero(cm.argument_size()),s); }
    else if (cm._terms.front().index().degree()==0) { cm._terms.front().coefficient()+=s; }
    else { cm._terms.prepend(MultiIndex::zero(cm.argument_size()),s); }
    return cm;
}

template<class X> auto
MultivariateChebyshevPolynomial<X>::apply(Mul, ChebyshevPolynomial<X> cm, Scalar<X> const& s) -> ChebyshevPolynomial<X> {
    for(auto term : cm._terms) { term.coefficient() *= s; }
    return cm;
}

template<class X> auto
MultivariateChebyshevPolynomial<X>::apply(Add, ChebyshevPolynomial<X> const& cm1, ChebyshevPolynomial<X> const& cm2) -> ChebyshevPolynomial<X> {
    ARIADNE_PRECONDITION(cm1.argument_size()==cm2.argument_size());
    ChebyshevPolynomial<X> cmr(cm1.argument_size(),min(cm1.precision(),cm2.precision()));
    auto iter1=cm1._terms.begin();
    auto iter2=cm2._terms.begin();
    while(iter1!=cm1._terms.end() && iter2!=cm2._terms.end()) {
        if(iter1->index()==iter2->index()) {
            cmr._terms.append(iter1->index(),iter1->coefficient()+iter2->coefficient());
            ++iter1; ++iter2;
        } else if(reverse_lexicographic_less(iter1->index(),iter2->index())) {
            cmr._terms.append(iter1->index(),iter1->coefficient());
            ++iter1;
        } else { // iter1->index()>iter2->index()
            cmr._terms.append(iter2->index(),iter2->coefficient());
            ++iter2;
        }
    }
    while(iter1!=cm1._terms.end()) {
        cmr._terms.append(iter1->index(),iter1->coefficient());
        ++iter1;
    }
    while(iter2!=cm2._terms.end()) {
        cmr._terms.append(iter2->index(),iter2->coefficient());
        ++iter2;
    }
    return cmr;
}

template<class X> auto
MultivariateChebyshevPolynomial<X>::apply(Sub, ChebyshevPolynomial<X> const& cm1, ChebyshevPolynomial<X> const& cm2) -> ChebyshevPolynomial<X> {
    return add(cm1,neg(cm2));
}

template<class X> auto
MultivariateChebyshevPolynomial<X>::apply(Mul, ChebyshevPolynomial<X> const& cm1, ChebyshevPolynomial<X> const& cm2) -> ChebyshevPolynomial<X> {
    ARIADNE_PRECONDITION(cm1.argument_size()==cm2.argument_size());
    SizeType const n=cm1.argument_size();
    PR pr=min(cm1.precision(),cm2.precision());
    ChebyshevPolynomial<X> result(n,pr);
    X const& z=result.zero_coefficient();

    // T_i(t) T_j(t) = (T_{i+j}(t) + T_{|i-j|}(t))/2.
    // If either degree is zero, the two contributions coincide.
    std::map<MultiIndex,X,ReverseLexicographicLess> coefficients;
    for(auto term1 : cm1._terms) {
        for(auto term2 : cm2._terms) {
            X product=term1.coefficient()*term2.coefficient();
            if(same(product,z)) { continue; }
            auto a=term1.index();
            auto b=term2.index();
            MultiIndex index(n);
            auto accumulate = [&](auto&& self, SizeType k, X const& c) -> Void {
                if(k==n) {
                    auto [iter,inserted]=coefficients.try_emplace(index,c);
                    if(not inserted) { iter->second+=c; }
                    return;
                }
                DegreeType i=a[k];
                DegreeType j=b[k];
                if(i==0u || j==0u) {
                    index[k]=static_cast<DegreeType>(i+j);
                    self(self,k+1u,c);
                } else {
                    X half=c/2;
                    index[k]=static_cast<DegreeType>(i+j);
                    self(self,k+1u,half);
                    index[k]=static_cast<DegreeType>(i>=j ? i-j : j-i);
                    self(self,k+1u,half);
                }
            };
            accumulate(accumulate,0u,product);
        }
    }

    for(auto const& [index,coefficient] : coefficients) {
        if(not same(coefficient,z)) { result._terms.append(index,coefficient); }
    }
    return result;
}


template<class X> auto
MultivariateChebyshevPolynomial<X>::apply(Sqr, ChebyshevPolynomial<X> cm) -> ChebyshevPolynomial<X> {
    return cm*cm;
}

template<class X> struct PartialChebyshevPolynomial {
    typedef typename MultivariateChebyshevPolynomial<X>::ConstIterator ConstIterator;
    ConstIterator _from, _to; SizeType _k;
    PartialChebyshevPolynomial(ConstIterator from, SizeType k, ConstIterator end)
        : _from(from), _to(from), _k(k) { while(_to!=end && equal_up_to(_from->index(),_to->index(),k)) { ++_to; } }
};

template<class X> class ChebyshevValues {
    mutable List<X> _lst;
  public:
    ChebyshevValues(X const& x) : _lst() { _lst.reserve(16); _lst.append(X(1,x.precision())); _lst.append(x); }
    X const& operator[] (SizeType i) const {
        for(SizeType j=_lst.size(); j<=i; ++j) { _lst.append(2*_lst[1]*_lst[j-1]-_lst[j-2]); } return _lst[i]; }
};

inline Boolean equal_up_to(MultiIndex const& a1, MultiIndex const& a2, SizeType k) {
    for (SizeType i=0; i!=k; ++i) { if (a1[i]!=a2[i]) { return false; } } return true;
}

template<class X> auto
MultivariateChebyshevPolynomial<X>::operator() (Vector<Scalar<X>> const& x) const -> Scalar<X> {
    ARIADNE_PRECONDITION(x.size()==this->argument_size());

    X r(x.zero_element());
    if(this->_terms.empty()) { return r; }

    List<ChebyshevValues<X>> values;
    values.reserve(x.size());
    for(SizeType i=0u; i!=x.size(); ++i) {
        values.append(ChebyshevValues<X>(x[i]));
    }

    for(SizeType j=0u; j!=this->_terms.size(); ++j) {
        auto const a=this->_terms._indices[j];
        X term=this->_terms._coefficients[j];
        for(SizeType i=0u; i!=x.size(); ++i) {
            term*=values[i][a[i]];
        }
        r+=term;
    }
    return r;
}



template<class X> auto
MultivariateChebyshevPolynomial<X>::_write(OutputStream& os) const -> OutputStream& {
    for (auto term_iter = this->_terms.begin(); term_iter!=this->_terms.end(); ++term_iter) {
        auto const& term = *term_iter;
        String cs = to_str(term.coefficient());
        if (cs.find('-')>=cs.find_first_of("0123456789") && term_iter!=this->_terms.begin()) { os << "+"; }
        os << cs;
        for (SizeType i=0; i!=this->argument_size(); ++i) { if (term.index()[i]!=0) { os << "*T" << term.index()[i] << "(x" << i << ")"; } }
    }
    return os;
}


} //namespace Ariadne



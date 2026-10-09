/***************************************************************************
 *            test_algebra.cpp
 *
 *  Copyright  2023  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */


#include "algebra/algebra.hpp"
#include "algebra/differential.hpp"

#include "utility/test.hpp"

using namespace Ariadne;

namespace {

template<class X> class AlgebraInterfaceStub : public AlgebraInterface<X> {
  public:
    AlgebraInterface<X>* _copy() const override { return new AlgebraInterfaceStub(*this); }
    AlgebraInterface<X>* _create_copy() const override { return new AlgebraInterfaceStub(*this); }
    AlgebraInterface<X>* _create_zero() const override { return new AlgebraInterfaceStub(*this); }
    AlgebraInterface<X>* _create_constant(X const&) const override { return new AlgebraInterfaceStub(*this); }

    Void _iadd(X const&) override { }
    Void _imul(X const&) override { }
    Void _isma(X const&, AlgebraInterface<X> const&) override { }
    Void _ifma(AlgebraInterface<X> const&, AlgebraInterface<X> const&) override { }

    AlgebraInterface<X>* _apply(UnaryRingOperator) const override { return new AlgebraInterfaceStub(*this); }
    AlgebraInterface<X>* _apply(BinaryRingOperator,AlgebraInterface<X> const&) const override { return new AlgebraInterfaceStub(*this); }
    AlgebraInterface<X>* _apply(BinaryFieldOperator,X const&) const override { return new AlgebraInterfaceStub(*this); }
    AlgebraInterface<X>* _rapply(BinaryRingOperator,X const&) const override { return new AlgebraInterfaceStub(*this); }
    AlgebraInterface<X>* _apply(GradedRingOperator,Nat) const override { return new AlgebraInterfaceStub(*this); }

    OutputStream& _write(OutputStream& os) const override { return os << "AlgebraInterfaceStub"; }
};

using ValidatedAlgebraInterfaceStub=AlgebraInterfaceStub<ValidatedNumber>;

template<class X> class TranscendentalAlgebraInterfaceStub : public TranscendentalAlgebraInterface<X> {
  public:
    TranscendentalAlgebraInterface<X>* _create_copy() const override { return new TranscendentalAlgebraInterfaceStub(*this); }
    TranscendentalAlgebraInterface<X>* _create_zero() const override { return new TranscendentalAlgebraInterfaceStub(*this); }
    TranscendentalAlgebraInterface<X>* _create_constant(X const&) const override { return new TranscendentalAlgebraInterfaceStub(*this); }

    TranscendentalAlgebraInterface<X>* _apply(BinaryFieldOperator,TranscendentalAlgebraInterface<X> const&) const override { return new TranscendentalAlgebraInterfaceStub(*this); }
    TranscendentalAlgebraInterface<X>* _apply(BinaryFieldOperator,X const&) const override { return new TranscendentalAlgebraInterfaceStub(*this); }
    TranscendentalAlgebraInterface<X>* _rapply(BinaryFieldOperator,X const&) const override { return new TranscendentalAlgebraInterfaceStub(*this); }
    TranscendentalAlgebraInterface<X>* _apply(UnaryTranscendentalOperator) const override { return new TranscendentalAlgebraInterfaceStub(*this); }
    TranscendentalAlgebraInterface<X>* _apply(GradedFieldOperator,Int) const override { return new TranscendentalAlgebraInterfaceStub(*this); }

    OutputStream& _write(OutputStream& os) const override { return os << "TranscendentalAlgebraInterfaceStub"; }
};

} // namespace

class TestAlgebra {
    using X=FloatDPApproximation;
    using DX=Differential<X>;
    using AX=Algebra<X>;
    using TAX=TranscendentalAlgebra<X>;
    using EAX=ElementaryAlgebra<X>;
  public:
    TestAlgebra() { }

    void test() const {
        DP pr;
        SizeType as=3;
        SizeType ind=1;
        DegreeType deg=4;
        ARIADNE_TEST_CONSTRUCT(X,cx,(3.0_x,pr));
        ARIADNE_TEST_NAMED_CONSTRUCT(DX,dx,variable(as,deg,cx,ind));
        ARIADNE_TEST_CONSTRUCT(AX,ax,(dx));
        AX rx=ax;
        ARIADNE_TEST_ASSIGN(ax,cx);
        ax=cx;
        ARIADNE_TEST_EXECUTE(ax+ax);
        ax+ax;
        ax-ax;
        ax+=ax;
        rx-=ax;
        ax+=cx;
        ax-=cx;
        ax=-ax;
        ax=ax+ax;
        ax=ax*ax;
        ax=cx+ax;
        ax=cx*ax;
        ax=ax+cx;
        ax=ax-cx;
        ax=ax*cx;
        ARIADNE_TEST_EXECUTE(ax=ax/cx);
        ARIADNE_TEST_ASSIGN(dx,ax.extract<DX>());
        auto algebra_wrapper=ax.extract<AlgebraWrapper<DX,X>>();
        ARIADNE_TEST_PRINT(algebra_wrapper);
        Algebra<X> interface_algebra(new AlgebraInterfaceStub<X>());
        ARIADNE_TEST_EXECUTE(interface_algebra.extract<AlgebraInterfaceStub<X>>());
        ARIADNE_TEST_FAIL(ax.extract<AlgebraInterfaceStub<X>>());
        ARIADNE_TEST_FAIL(interface_algebra.extract<DX>());

        Algebra<ValidatedNumber> validated_algebra(new ValidatedAlgebraInterfaceStub());
        AlgebraInterface<ValidatedNumber> const& validated_interface=validated_algebra;
        (void)validated_interface;
        ARIADNE_TEST_EXECUTE(validated_algebra.create());
        ARIADNE_TEST_EXECUTE(validated_algebra.clone());
        ARIADNE_TEST_EXECUTE(validated_algebra.create_zero());
        auto validated_body=validated_algebra.extract<ValidatedAlgebraInterfaceStub>();
        ARIADNE_TEST_PRINT(validated_body);

        ARIADNE_TEST_CONSTRUCT(TAX,tax,(dx));
        auto transcendental_zero=tax.reference()._create_zero();
        delete transcendental_zero;
        auto transcendental_constant=tax.reference()._create_constant(cx);
        delete transcendental_constant;
        auto transcendental_copy=tax.reference()._create_copy();
        delete transcendental_copy;
        auto transcendental_rdiv=tax.reference()._rapply(BinaryFieldOperator(Div()),cx);
        delete transcendental_rdiv;
        ARIADNE_TEST_ASSIGN(tax,div(tax,tax+cx));
        ARIADNE_TEST_ASSIGN(tax,rec(tax));

        ARIADNE_TEST_EXECUTE(hlf(tax));
        ARIADNE_TEST_EXECUTE(div(cx,tax));
        ARIADNE_TEST_EXECUTE(pow(tax,Int(2)));
        ARIADNE_TEST_EXECUTE(pow(tax,Int(-1)));
        ARIADNE_TEST_EXECUTE(sqrt(tax));
        ARIADNE_TEST_EXECUTE(log(tax));
        ARIADNE_TEST_EXECUTE(sin(tax));
        ARIADNE_TEST_EXECUTE(cos(tax));
        ARIADNE_TEST_EXECUTE(tan(tax));
        ARIADNE_TEST_EXECUTE(tanh(tax));
        ARIADNE_TEST_EXECUTE(atan(tax));

        ARIADNE_TEST_EXECUTE(GradedAlgebraOperations<DX>::apply(Hlf(),dx));
        ARIADNE_TEST_EXECUTE(GradedAlgebraOperations<DX>::apply(Sqr(),dx));
        ARIADNE_TEST_EXECUTE(GradedAlgebraOperations<DX>::apply(Div(),dx,dx));
        ARIADNE_TEST_EXECUTE(GradedAlgebraOperations<DX>::apply(Pow(),dx,Int(2)));
        ARIADNE_TEST_EXECUTE(GradedAlgebraOperations<DX>::apply(Pow(),dx,Int(-1)));
        ARIADNE_TEST_EXECUTE(GradedAlgebraOperations<DX>::apply(Exp(),dx));

        Factorial factorial(3u);
        FloatDPBounds factorial_value=factorial;
        ARIADNE_TEST_PRINT(factorial_value);
        ARIADNE_TEST_EXECUTE(rec(factorial));

        ARIADNE_TEST_ASSIGN(dx,tax.extract<DX>());
        auto transcendental_wrapper=tax.extract<TranscendentalAlgebraWrapper<DX,X>>();
        ARIADNE_TEST_PRINT(transcendental_wrapper);
        TranscendentalAlgebra<X> interface_tax(new TranscendentalAlgebraInterfaceStub<X>());
        auto copied_transcendental_interface=interface_tax.reference()._copy();
        delete copied_transcendental_interface;
        ARIADNE_TEST_EXECUTE(interface_tax.extract<TranscendentalAlgebraInterfaceStub<X>>());
        ARIADNE_TEST_FAIL(tax.extract<TranscendentalAlgebraInterfaceStub<X>>());
        ARIADNE_TEST_FAIL(interface_tax.extract<DX>());

    }
};

Int main() {
    TestAlgebra().test();
    return ARIADNE_TEST_FAILURES;
}

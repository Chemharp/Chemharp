// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#ifndef CHEMFILES_SELECTION_EXPR_HPP
#define CHEMFILES_SELECTION_EXPR_HPP

#include <cassert>
#include <cstddef>
#include <cstdint>

#include <string>
#include <utility>
#include <vector>
#include <memory>
#include <functional>

#include "chemfiles/external/optional.hpp"
#include "chemfiles/selections/NumericValues.hpp"

namespace chemfiles {
class Frame;
class Match;
class Selection;

namespace selections {

using Variable = uint8_t;

/// Abstract base class for selectors in the selection AST
class Selector {
public:
    /// Pretty-printing of this selector. The output should use a shift
    /// of `delta` spaces in case of multilines output.
    virtual std::string print(unsigned delta = 0) const = 0;
    /// Check if the `match` is valid in the given `frame`.
    virtual bool is_match(const Frame& frame, const Match& match) const = 0;
    /// Clear any cached data. This must be called before using the selection
    /// with a new frame
    virtual void clear() = 0;
    /// Optimize the AST corresponding to this Selector. Currently, this only
    /// perform constant propgations in mathematical expressions.
    virtual void optimize() {}

    Selector() = default;
    virtual ~Selector() = default;

    // Delete all move and copy constructors.
    // Selector should only be used behind an std::unique_ptr
    Selector(Selector&&) = delete;
    Selector& operator=(Selector&&) = delete;
    Selector(const Selector&) = delete;
    Selector& operator=(const Selector&) = delete;
};

using Ast = std::unique_ptr<Selector>;

/// Combine selections by using a logical `and` operation
class And final: public Selector {
public:
    And(Ast lhs, Ast rhs): lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}
    std::string print(unsigned delta) const override;
    bool is_match(const Frame& frame, const Match& match) const override;
    void clear() override;
private:
    Ast lhs_;
    Ast rhs_;
};

/// Combine selections by using a logical `or` operation
class Or final: public Selector {
public:
    Or(Ast lhs, Ast rhs): lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}
    std::string print(unsigned delta) const override;
    bool is_match(const Frame& frame, const Match& match) const override;
    void clear() override;
private:
    Ast lhs_;
    Ast rhs_;
};

/// Unary negation of a selection
class Not final: public Selector {
public:
    explicit Not(Ast ast): ast_(std::move(ast)) {}
    std::string print(unsigned delta) const override;
    bool is_match(const Frame& frame, const Match& match) const override;
    void clear() override;
private:
    Ast ast_;
};

/// Selection matching all atoms
class All final: public Selector {
public:
    All() = default;
    std::string print(unsigned delta) const override;
    bool is_match(const Frame& frame, const Match& match) const override;
    void clear() override {}
};

/// Selection matching no atoms
class None final: public Selector {
public:
    None() = default;
    std::string print(unsigned delta) const override;
    bool is_match(const Frame& frame, const Match& match) const override;
    void clear() override {}
};

/// Selection based on boolean properties
class BoolProperty final: public Selector {
public:
    BoolProperty(std::string property, Variable argument):
        property_(std::move(property)), argument_(argument) {}
    std::string print(unsigned delta) const override;
    bool is_match(const Frame& frame, const Match& match) const override;
    void clear() override {}

private:
    std::string property_;
    /// Which atom in the candidate match are we checking?
    Variable argument_;
};

/// A sub-selection for use in boolean selectors
class SubSelection {
public:
    /// Create a sub-selection from a variable
    SubSelection(Variable variable);
    /// Create a sub-selection from an AST
    SubSelection(std::string selection);

    SubSelection(SubSelection&&) noexcept;
    SubSelection& operator=(SubSelection&&) noexcept;
    ~SubSelection();

    /// Evaluate the sub-selection and return the list of matching atoms
    const std::vector<size_t>& eval(const Frame& frame, const Match& match) const;
    /// Pretty-print the sub-selection
    std::string print() const;
    /// Clear cached data
    void clear();

    bool is_variable() const {
        return selection_ == nullptr;
    }

private:
    /// Possible selection. If this is nullptr, then the variable_ is set.
    std::unique_ptr<Selection> selection_;
    /// Variable to use if selection_ is nullptr
    Variable variable_;
    /// Cache matches for the selection on the first call to eval
    mutable std::vector<size_t> matches_;
    /// Did we update the cached matches ?
    mutable bool updated_ = false;
};

/// Checking if two atoms are bonded together
class IsBonded final: public Selector {
public:
    IsBonded(SubSelection i, SubSelection j): i_(std::move(i)), j_(std::move(j)) {}
    std::string print(unsigned delta) const override;
    bool is_match(const Frame& frame, const Match& match) const override;
    void clear() override;
private:
    SubSelection i_;
    SubSelection j_;
};

/// Checking if three atoms are bonded together to form an angle
class IsAngle final: public Selector {
public:
    IsAngle(SubSelection i, SubSelection j, SubSelection k):
        i_(std::move(i)), j_(std::move(j)), k_(std::move(k)) {}
    std::string print(unsigned delta) const override;
    bool is_match(const Frame& frame, const Match& match) const override;
    void clear() override;
private:
    SubSelection i_;
    SubSelection j_;
    SubSelection k_;
};

/// Checking if four atoms are bonded together to form a dihedral angle
class IsDihedral final: public Selector {
public:
    IsDihedral(SubSelection i, SubSelection j, SubSelection k, SubSelection m):
        i_(std::move(i)), j_(std::move(j)), k_(std::move(k)), m_(std::move(m)) {}
    std::string print(unsigned delta) const override;
    bool is_match(const Frame& frame, const Match& match) const override;
    void clear() override;
private:
    SubSelection i_;
    SubSelection j_;
    SubSelection k_;
    SubSelection m_;
};

/// Checking if four atoms are bonded together to form an improper dihedral angle
class IsImproper final: public Selector {
public:
    IsImproper(SubSelection i, SubSelection j, SubSelection k, SubSelection m):
        i_(std::move(i)), j_(std::move(j)), k_(std::move(k)), m_(std::move(m)) {}
    std::string print(unsigned delta) const override;
    bool is_match(const Frame& frame, const Match& match) const override;
    void clear() override;
private:
    SubSelection i_;
    SubSelection j_;
    SubSelection k_;
    SubSelection m_;
};

/// Abstract base class for string selector
class StringSelector: public Selector {
public:
    StringSelector(std::string value, bool equals, Variable argument):
        value_(std::move(value)), equals_(equals), argument_(argument)
    {
        assert(argument <= 3 && "argument must be less than 3 in SingleSelector");
    }
    ~StringSelector() override = default;

    /// Get the value for the atom at index `i` in the `frame`
    virtual const std::string& value(const Frame& frame, size_t i) const = 0;
    /// Get the property name
    virtual std::string name() const = 0;

    bool is_match(const Frame& frame, const Match& match) const final;
    std::string print(unsigned delta) const final;

private:
    /// The value to check against
    std::string value_;
    /// Are we checking for equality or inequality?
    bool equals_;
    /// Which atom in the candidate match are we checking?
    Variable argument_;
};

/// Selection based on string properties
class StringProperty final: public StringSelector {
public:
    StringProperty(std::string property, std::string value, bool equals, Variable argument):
        StringSelector(std::move(value), equals, argument),
        property_(std::move(property)) {}

    const std::string& value(const Frame& frame, size_t i) const override;
    std::string name() const override;
    void clear() override {}

private:
    std::string property_;
};

/// Select atoms using their type
class Type final: public StringSelector {
public:
    Type(std::string value, bool equals, Variable argument):
        StringSelector(std::move(value), equals, argument) {}

    std::string name() const override;
    const std::string& value(const Frame& frame, size_t i) const override;
    void clear() override {}
};

/// Select atoms using their name
class Name final: public StringSelector {
public:
    Name(std::string value, bool equals, Variable argument):
        StringSelector(std::move(value), equals, argument) {}

    std::string name() const override;
    const std::string& value(const Frame& frame, size_t i) const override;
    void clear() override {}
};

/// Select atoms using their residue name
class Resname final: public StringSelector {
public:
    Resname(std::string value, bool equals, Variable argument):
        StringSelector(std::move(value), equals, argument) {}

    std::string name() const override;
    const std::string& value(const Frame& frame, size_t i) const override;
    void clear() override {}
};

class MathExpr;
using MathAst = std::unique_ptr<MathExpr>;

/// Expression for math selectors.
///
/// Math selectors look like `<lhs> <op> <rhs>`, where `<lhs>` and `<rhs>`
/// evaluate to numeric values, and `<op>` can be any comparison operator (==,
/// <, >=, !=, ...). This makes the whole selector evaluate to a boolean value.
///
/// It is possible for either or both `<lhs>` and `<rhs>` to evaluate to
/// multiple numeric values (for example `distance(#1, name O)`), in which case
/// the math selector with evaluate to `true` if ANY of the comparison evaluate
/// to `true`.
class Math final: public Selector {
public:
    enum class Operator {
        EQUAL,
        NOT_EQUAL,
        LESS,
        LESS_EQUAL,
        GREATER,
        GREATER_EQUAL,
    };

    Math(Operator op, MathAst lhs, MathAst rhs): op_(op), lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}

    bool is_match(const Frame& frame, const Match& match) const override;
    void optimize() override;
    std::string print(unsigned delta) const override;
    void clear() override;

private:
    Operator op_;
    MathAst lhs_;
    MathAst rhs_;
};

/// Abstract base class for mathematical expressions
class MathExpr {
public:
    MathExpr() = default;
    virtual ~MathExpr() = default;

    // Delete all move and copy constructors.
    // MathExpr should only be used behind an std::unique_ptr
    MathExpr(MathExpr&&) = delete;
    MathExpr& operator=(MathExpr&&) = delete;
    MathExpr(const MathExpr&) = delete;
    MathExpr& operator=(const MathExpr&) = delete;

    /// Evaluate the expression and get the all the matching values
    virtual NumericValues eval(const Frame& frame, const Match& match) const = 0;

    /// Propagate all constants in this sub ast, and return the corresponding
    /// value if possible.
    ///
    /// This currently only support single values, i.e. optimizing `3 + 4` to 7.
    virtual optional<double> optimize() = 0;

    /// Clear any cached data
    virtual void clear() = 0;

    /// Pretty-print the expression
    virtual std::string print() const = 0;
};

// Addition
class Add final: public MathExpr {
public:
    Add(MathAst lhs, MathAst rhs): lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override;
    std::string print() const override;
    void clear() override;

private:
    MathAst lhs_;
    MathAst rhs_;
};

// Substraction
class Sub final: public MathExpr {
public:
    Sub(MathAst lhs, MathAst rhs): lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override;
    std::string print() const override;
    void clear() override;

private:
    MathAst lhs_;
    MathAst rhs_;
};

// Multiplication
class Mul final: public MathExpr {
public:
    Mul(MathAst lhs, MathAst rhs): lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override;
    std::string print() const override;
    void clear() override;

private:
    MathAst lhs_;
    MathAst rhs_;
};

// Division
class Div final: public MathExpr {
public:
    Div(MathAst lhs, MathAst rhs): lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override;
    std::string print() const override;
    void clear() override;

private:
    MathAst lhs_;
    MathAst rhs_;
};

/// Power raising
class Pow final: public MathExpr {
public:
    Pow(MathAst lhs, MathAst rhs): lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override;
    std::string print() const override;
    void clear() override;

private:
    MathAst lhs_;
    MathAst rhs_;
};

/// Unary minus operator
class Neg final: public MathExpr {
public:
    Neg(MathAst ast): ast_(std::move(ast)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override;
    std::string print() const override;
    void clear() override;

private:
    MathAst ast_;
};

/// Modulo (remainder of Euclidean division) operation
class Mod final: public MathExpr {
public:
    Mod(MathAst lhs, MathAst rhs): lhs_(std::move(lhs)), rhs_(std::move(rhs)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override;
    std::string print() const override;
    void clear() override;

private:
    MathAst lhs_;
    MathAst rhs_;
};

/// Function call
class Function final: public MathExpr {
public:
    Function(std::function<double(double)> fn, std::string name, MathAst ast):
        fn_(std::move(fn)), name_(std::move(name)), ast_(std::move(ast)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override;
    std::string print() const override;
    void clear() override;

private:
    std::function<double(double)> fn_;
    std::string name_;
    MathAst ast_;
};

/// Statically known number (either literal number or the result of
/// optimization of literal numbers).
class Number final: public MathExpr {
public:
    Number(double value): value_(value) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override;
    std::string print() const override;
    void clear() override {}

private:
    double value_;
};

/// Compute the distance between atoms
class Distance final: public MathExpr {
public:
    Distance(SubSelection i, SubSelection j): i_(std::move(i)), j_(std::move(j)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override {
        return nullopt;
    }
    std::string print() const override;
    void clear() override {}

private:
    SubSelection i_;
    SubSelection j_;
};

/// Compute the angle between three atoms
class Angle final: public MathExpr {
public:
    Angle(SubSelection i, SubSelection j, SubSelection k):
        i_(std::move(i)), j_(std::move(j)), k_(std::move(k)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override {
        return nullopt;
    }
    std::string print() const override;
    void clear() override {}

private:
    SubSelection i_;
    SubSelection j_;
    SubSelection k_;
};

/// Compute the dihedral angle between four atoms
class Dihedral final: public MathExpr {
public:
    Dihedral(SubSelection i, SubSelection j, SubSelection k, SubSelection m):
        i_(std::move(i)), j_(std::move(j)), k_(std::move(k)), m_(std::move(m)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override {
        return nullopt;
    }
    std::string print() const override;
    void clear() override {}

private:
    SubSelection i_;
    SubSelection j_;
    SubSelection k_;
    SubSelection m_;
};

/// Compute the out of plane distance between four atoms
class OutOfPlane final: public MathExpr {
public:
    OutOfPlane(SubSelection i, SubSelection j, SubSelection k, SubSelection m):
        i_(std::move(i)), j_(std::move(j)), k_(std::move(k)), m_(std::move(m)) {}

    NumericValues eval(const Frame& frame, const Match& match) const override;
    optional<double> optimize() override {
        return nullopt;
    }
    std::string print() const override;
    void clear() override {}

private:
    SubSelection i_;
    SubSelection j_;
    SubSelection k_;
    SubSelection m_;
};

/// Abstract base class for numeric selectors
class NumericSelector: public MathExpr {
public:
    NumericSelector(Variable argument): argument_(argument) {}
    ~NumericSelector() override = default;

    NumericValues eval(const Frame& frame, const Match& match) const final;
    optional<double> optimize() final;
    std::string print() const final;

    /// Get the value for the atom at index `i` in the `frame`
    virtual double value(const Frame& frame, size_t i) const = 0;
    /// Get the name of the selector
    virtual std::string name() const = 0;

private:
    /// Which atom in the candidate match are we checking?
    Variable argument_;
};

/// Select atoms using a given double property in the frame.
class NumericProperty final: public NumericSelector {
public:
    NumericProperty(std::string property, Variable argument): NumericSelector(argument), property_(std::move(property)) {}
    std::string name() const override;
    double value(const Frame& frame, size_t i) const override;
    void clear() override {}

private:
    std::string property_;
};


/// Select atoms using their index in the frame.
class Index final: public NumericSelector {
public:
    Index(Variable argument): NumericSelector(argument) {}
    std::string name() const override;
    double value(const Frame& frame, size_t i) const override;
    void clear() override {}
};

/// Select atoms using their residue id (residue number)
class Resid final: public NumericSelector {
public:
    Resid(Variable argument): NumericSelector(argument) {}
    std::string name() const override;
    double value(const Frame& frame, size_t i) const override;
    void clear() override {}
};

/// Select atoms using their mass.
class Mass final: public NumericSelector {
public:
    Mass(Variable argument): NumericSelector(argument) {}
    std::string name() const override;
    double value(const Frame& frame, size_t i) const override;
    void clear() override {}
};

enum class Coordinate {
    X = 0,
    Y = 1,
    Z = 2,
};

/// Select atoms using their position in space. The selection can be created by
/// `x <op> <val>`, `y <op> <val>` or `z <op> <val>`, depending on the component
/// of the position to use.
class Position final: public NumericSelector {
public:
    Position(Variable argument, Coordinate coordinate): NumericSelector(argument), coordinate_(coordinate) {}
    std::string name() const override;
    double value(const Frame& frame, size_t i) const override;
    void clear() override {}

private:
    Coordinate coordinate_;
};

/// Select atoms using their velocity. The selection can be created by `vx <op>
/// <val>`, `vy <op> <val>` or `vz <op> <val>`, depending on the component of
/// the velocity to use.
class Velocity final: public NumericSelector {
public:
    Velocity(Variable argument, Coordinate coordinate): NumericSelector(argument), coordinate_(coordinate) {}
    std::string name() const override;
    double value(const Frame& frame, size_t i) const override;
    void clear() override {}

private:
    Coordinate coordinate_;
};

}} // namespace chemfiles && namespace selections

#endif

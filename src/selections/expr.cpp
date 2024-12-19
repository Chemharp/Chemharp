// Chemfiles, a modern library for chemistry file reading and writing
// Copyright (C) Guillaume Fraux and contributors -- BSD license

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cassert>

#include <array>
#include <string>
#include <utility>
#include <vector>
#include <memory>
#include <algorithm>
#include <functional>

#include <fmt/core.h>
#include <fmt/format.h>

#include "chemfiles/error_fmt.hpp"
#include "chemfiles/unreachable.hpp"
#include "chemfiles/external/optional.hpp"

#include "chemfiles/Atom.hpp"
#include "chemfiles/Frame.hpp"
#include "chemfiles/Property.hpp"
#include "chemfiles/Residue.hpp"
#include "chemfiles/Selection.hpp"
#include "chemfiles/Topology.hpp"
#include "chemfiles/Connectivity.hpp"

#include "chemfiles/selections/expr.hpp"
#include "chemfiles/selections/lexer.hpp"
#include "chemfiles/selections/NumericValues.hpp"

using namespace chemfiles;
using namespace chemfiles::selections;

static const std::string EMPTY_STRING;

static std::string kind_as_string(Property::Kind kind) {
    switch (kind) {
    case Property::BOOL:
        return "bool";
    case Property::DOUBLE:
        return "double";
    case Property::STRING:
        return "string";
    case Property::VECTOR3D:
        return "Vector3D";
    default:
        unreachable();
    }
}

SubSelection::SubSelection(SubSelection&&) noexcept = default;
SubSelection& SubSelection::operator=(SubSelection&&) noexcept = default;
SubSelection::~SubSelection() = default;

SubSelection::SubSelection(Variable variable): selection_(nullptr), variable_(variable) {
    matches_.resize(1);
}

SubSelection::SubSelection(std::string selection):
    selection_(std::make_unique<Selection>(std::move(selection))), variable_(UINT8_MAX)
{
    assert(selection_->size() == 1);
}

const std::vector<size_t>& SubSelection::eval(const Frame& frame, const Match& match) const {
    if (is_variable()) {
        matches_[0] = match[variable_];
    } else {
        if (!updated_) {
            matches_ = selection_->list(frame);
            updated_ = true;
        }
    }
    return matches_;
}

std::string SubSelection::print() const {
    if (is_variable()) {
        return fmt::format("#{}", variable_ + 1);
    } else {
        return selection_->string();
    }
}

void SubSelection::clear() {
    if (!is_variable()) {
        matches_.clear();
        updated_ = false;
    }
}

std::string And::print(unsigned delta) const {
    auto lhs = lhs_->print(7);
    auto rhs = rhs_->print(7);
    return fmt::format("and -> {}\n{:{}}    -> {}", lhs, "", delta, rhs);
}

bool And::is_match(const Frame& frame, const Match& match) const {
    return lhs_->is_match(frame, match) && rhs_->is_match(frame, match);
}

void And::clear() {
    lhs_->clear();
    rhs_->clear();
}

std::string Or::print(unsigned delta) const {
    auto lhs = lhs_->print(6);
    auto rhs = rhs_->print(6);
    return fmt::format("or -> {}\n{:{}}   -> {}", lhs, "", delta, rhs);
}

bool Or::is_match(const Frame& frame, const Match& match) const {
    return lhs_->is_match(frame, match) || rhs_->is_match(frame, match);
}

void Or::clear() {
    lhs_->clear();
    rhs_->clear();
}

std::string Not::print(unsigned /*unused*/) const {
    return "not " + ast_->print(4);
}

bool Not::is_match(const Frame& frame, const Match& match) const {
    return !ast_->is_match(frame, match);
}

void Not::clear() {
    ast_->clear();
}

std::string All::print(unsigned /*unused*/) const {
    return "all";
}

bool All::is_match(const Frame& /*unused*/, const Match& /*unused*/) const {
    return true;
}

std::string None::print(unsigned /*unused*/) const {
    return "none";
}

bool None::is_match(const Frame& /*unused*/, const Match& /*unused*/) const {
    return false;
}

std::string BoolProperty::print(unsigned /*unused*/) const {
    if (is_ident(property_)) {
        return fmt::format("[{}](#{})", property_, argument_ + 1);
    } else {
        return fmt::format("[\"{}\"](#{})", property_, argument_ + 1);
    }
}

bool BoolProperty::is_match(const Frame& frame, const Match& match) const {
    const auto& property = frame[match[argument_]].get(property_);
    if (property) {
        if (property->kind() == Property::BOOL) {
            return property->as_bool();
        } else {
            throw selection_error(
                "invalid type for property [{}] on atom {}: expected bool, got {}",
                property_, match[argument_], kind_as_string(property->kind())
            );
        }
    } else {
        auto residue = frame.topology().residue_for_atom(match[argument_]);
        if (residue) {
            const auto& resProperty = residue->get(property_);
            if(!resProperty) {
                return false;
            }
            if (resProperty->kind() == Property::BOOL) {
                return resProperty->as_bool();
            } else {
                throw selection_error(
                    "invalid type for property [{}] on the residue containing atom {}: expected bool, got {}",
                    property_, match[argument_], kind_as_string(resProperty->kind())
                );
            }
        }
        // No property with the given name
        return false;
    }
}

std::string IsBonded::print(unsigned /*unused*/) const {
    return fmt::format("is_bonded({}, {})", i_.print(), j_ .print());
}

bool IsBonded::is_match(const Frame& frame, const Match& match) const {
    const auto& bonds = frame.topology().bonds();
    for (auto i: i_.eval(frame, match)) {
        for (auto j: j_.eval(frame, match)) {
            if (i == j) {
                // This is not a valid bond
                continue;
            }

            auto candidate = chemfiles::Bond(i, j);
            auto it = std::lower_bound(bonds.begin(), bonds.end(), candidate);
            if (it != bonds.end() && *it == candidate) {
                return true;
            }
        }
    }
    return false;
}

void IsBonded::clear() {
    i_.clear();
    j_.clear();
}

std::string IsAngle::print(unsigned /*unused*/) const {
    return fmt::format("is_angle({}, {}, {})", i_.print(), j_.print(), k_.print());
}

bool IsAngle::is_match(const Frame& frame, const Match& match) const {
    const auto& angles = frame.topology().angles();
    for (auto i: i_.eval(frame, match)) {
        for (auto j: j_.eval(frame, match)) {
            for (auto k: k_.eval(frame, match)) {
                if (i == j || j == k || i == k) {
                    // This is not a valid angle
                    continue;
                }

                auto candidate = chemfiles::Angle(i, j, k);
                auto it = std::lower_bound(angles.begin(), angles.end(), candidate);
                if (it != angles.end() && *it == candidate) {
                    return true;
                }
            }
        }
    }
    return false;
}

void IsAngle::clear() {
    i_.clear();
    j_.clear();
    k_.clear();
}

std::string IsDihedral::print(unsigned /*unused*/) const {
    return fmt::format("is_dihedral({}, {}, {}, {})", i_.print(), j_.print(), k_.print(), m_.print());
}

bool IsDihedral::is_match(const Frame& frame, const Match& match) const {
    const auto& dihedrals = frame.topology().dihedrals();
    for (auto i: i_.eval(frame, match)) {
        for (auto j: j_.eval(frame, match)) {
            for (auto k: k_.eval(frame, match)) {
                for (auto m: m_.eval(frame, match)) {
                    if (i == j || j == k || k == m || i == k || j == m || i == m) {
                        // This is not a valid angle
                        continue;
                    }

                    auto candidate = chemfiles::Dihedral(i, j, k, m);
                    auto it = std::lower_bound(dihedrals.begin(), dihedrals.end(), candidate);
                    if (it != dihedrals.end() && *it == candidate) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

void IsDihedral::clear() {
    i_.clear();
    j_.clear();
    k_.clear();
    m_.clear();
}

std::string IsImproper::print(unsigned /*unused*/) const {
    return fmt::format("is_improper({}, {}, {}, {})", i_.print(), j_.print(), k_.print(), m_.print());
}

bool IsImproper::is_match(const Frame& frame, const Match& match) const {
    const auto& impropers = frame.topology().impropers();
    for (auto i: i_.eval(frame, match)) {
        for (auto j: j_.eval(frame, match)) {
            for (auto k: k_.eval(frame, match)) {
                for (auto m: m_.eval(frame, match)) {
                    if (i == j || j == k || k == m || i == k || j == m || i == m) {
                        // This is not a valid angle
                        continue;
                    }

                    auto candidate = chemfiles::Improper(i, j, k, m);
                    auto it = std::lower_bound(impropers.begin(), impropers.end(), candidate);
                    if (it != impropers.end() && *it == candidate) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

void IsImproper::clear() {
    i_.clear();
    j_.clear();
    k_.clear();
    m_.clear();
}

std::string StringSelector::print(unsigned /*unused*/) const {
    const auto *op = equals_ ? "==" : "!=";
    if (is_ident(value_)) {
        return fmt::format("{}(#{}) {} {}", name(), argument_ + 1, op, value_);
    } else {
        return fmt::format("{}(#{}) {} \"{}\"", name(), argument_ + 1, op, value_);
    }
}

bool StringSelector::is_match(const Frame& frame, const Match& match) const {
    return (this->value(frame, match[argument_]) == value_) == equals_;
}

std::string StringProperty::name() const {
    if (is_ident(property_)) {
        return "[" + property_ + "]";
    } else {
        return "[\"" + property_ + "\"]";
    }
}

const std::string& StringProperty::value(const Frame& frame, size_t i) const {
    const auto& property = frame[i].get(property_);
    if (property) {
        if (property->kind() == Property::STRING) {
            return property->as_string();
        } else {
            throw selection_error(
                "invalid type for property [{}] on atom {}: expected string, got {}",
                property_, i, kind_as_string(property->kind())
            );
        }
    } else {
        auto residue = frame.topology().residue_for_atom(i);
        if (residue) {
            const auto& resProperty = residue->get(property_);
            if(!resProperty) {
                return EMPTY_STRING;
            }
            if (resProperty->kind() == Property::STRING) {
                return resProperty->as_string();
            } else {
                throw selection_error(
                    "invalid type for property [{}] on the residue containing atom {}: expected string, got {}",
                    property_, i, kind_as_string(resProperty->kind())
                );
            }
        }
        // No property with the given name
        return EMPTY_STRING;
    }
}

std::string Type::name() const {
    return "type";
}

const std::string& Type::value(const Frame& frame, size_t i) const {
    return frame[i].type();
}

std::string Name::name() const {
    return "name";
}

const std::string& Name::value(const Frame& frame, size_t i) const {
    return frame[i].name();
}

std::string Resname::name() const {
    return "resname";
}

const std::string& Resname::value(const Frame& frame, size_t i) const {
    auto residue = frame.topology().residue_for_atom(i);
    if (residue) {
        return residue->name();
    } else {
        return EMPTY_STRING;
    }
}

bool Math::is_match(const Frame& frame, const Match& match) const {
    std::function<bool(double, double)> operation;
    switch (op_) {
    case Math::Operator::EQUAL:
        operation = [](double l, double r){ return l == r; };
        break;
    case Math::Operator::NOT_EQUAL:
        operation = [](double l, double r){ return l != r; };
        break;
    case Math::Operator::LESS:
        operation = [](double l, double r){ return l < r; };
        break;
    case Math::Operator::LESS_EQUAL:
        operation = [](double l, double r){ return l <= r; };
        break;
    case Math::Operator::GREATER:
        operation = [](double l, double r){ return l > r; };
        break;
    case Math::Operator::GREATER_EQUAL:
        operation = [](double l, double r){ return l >= r; };
        break;
    default:
        unreachable();
    }

    auto lhs = lhs_->eval(frame, match);
    auto rhs = rhs_->eval(frame, match);
    for (auto left: lhs) {
        for (auto right: rhs) {
            if (operation(left, right)) {
                return true;
            }
        }
    }
    return false;
}

std::string Math::print(unsigned /*unused*/) const {
    std::string op;
    switch (op_) {
    case Math::Operator::EQUAL:
        op = " == ";
        break;
    case Math::Operator::NOT_EQUAL:
        op = " != ";
        break;
    case Math::Operator::LESS:
        op = " < ";
        break;
    case Math::Operator::LESS_EQUAL:
        op = " <= ";
        break;
    case Math::Operator::GREATER:
        op = " > ";
        break;
    case Math::Operator::GREATER_EQUAL:
        op = " >= ";
        break;
    default:
        unreachable();
    }

    return lhs_->print() + op + rhs_->print();
}

void Math::optimize() {
    auto lhs_opt = lhs_->optimize();
    if (lhs_opt) {
        lhs_ = std::make_unique<Number>(lhs_opt.value());
    }

    auto rhs_opt = rhs_->optimize();
    if (rhs_opt) {
        rhs_ = std::make_unique<Number>(rhs_opt.value());
    }
}

void Math::clear() {
    lhs_->clear();
    rhs_->clear();
}

NumericValues Add::eval(const Frame& frame, const Match& match) const {
    auto lhs = lhs_->eval(frame, match);
    auto rhs = rhs_->eval(frame, match);

    auto result = NumericValues();
    result.reserve(lhs.size() * rhs.size());
    for (auto left: lhs) {
        for (auto right: rhs) {
            result.push_back(left + right);
        }
    }

    return result;
}

optional<double> Add::optimize() {
    auto lhs_opt = lhs_->optimize();
    auto rhs_opt = rhs_->optimize();
    if (lhs_opt && rhs_opt) {
        return lhs_opt.value() + rhs_opt.value();
    } else if (lhs_opt) {
        lhs_ = std::make_unique<Number>(lhs_opt.value());
    } else if (rhs_opt) {
        rhs_ = std::make_unique<Number>(rhs_opt.value());
    }
    return nullopt;
}

std::string Add::print() const {
    return fmt::format("({} + {})", lhs_->print(), rhs_->print());
}

void Add::clear() {
    lhs_->clear();
    rhs_->clear();
}

NumericValues Sub::eval(const Frame& frame, const Match& match) const {
    auto lhs = lhs_->eval(frame, match);
    auto rhs = rhs_->eval(frame, match);

    auto result = NumericValues();
    result.reserve(lhs.size() * rhs.size());
    for (auto left: lhs) {
        for (auto right: rhs) {
            result.push_back(left - right);
        }
    }

    return result;
}

optional<double> Sub::optimize() {
    auto lhs_opt = lhs_->optimize();
    auto rhs_opt = rhs_->optimize();
    if (lhs_opt && rhs_opt) {
        return lhs_opt.value() - rhs_opt.value();
    } else if (lhs_opt) {
        lhs_ = std::make_unique<Number>(lhs_opt.value());
    } else if (rhs_opt) {
        rhs_ = std::make_unique<Number>(rhs_opt.value());
    }
    return nullopt;
}

std::string Sub::print() const {
    return fmt::format("({} - {})", lhs_->print(), rhs_->print());
}

void Sub::clear() {
    lhs_->clear();
    rhs_->clear();
}

NumericValues Mul::eval(const Frame& frame, const Match& match) const {
    auto lhs = lhs_->eval(frame, match);
    auto rhs = rhs_->eval(frame, match);

    auto result = NumericValues();
    result.reserve(lhs.size() * rhs.size());
    for (auto left: lhs) {
        for (auto right: rhs) {
            result.push_back(left * right);
        }
    }

    return result;
}

optional<double> Mul::optimize() {
    auto lhs_opt = lhs_->optimize();
    auto rhs_opt = rhs_->optimize();
    if (lhs_opt && rhs_opt) {
        return lhs_opt.value() * rhs_opt.value();
    } else if (lhs_opt) {
        lhs_ = std::make_unique<Number>(lhs_opt.value());
    } else if (rhs_opt) {
        rhs_ = std::make_unique<Number>(rhs_opt.value());
    }
    return nullopt;
}

std::string Mul::print() const {
    return fmt::format("({} * {})", lhs_->print(), rhs_->print());
}

void Mul::clear() {
    lhs_->clear();
    rhs_->clear();
}

NumericValues Div::eval(const Frame& frame, const Match& match) const {
    auto lhs = lhs_->eval(frame, match);
    auto rhs = rhs_->eval(frame, match);

    auto result = NumericValues();
    result.reserve(lhs.size() * rhs.size());
    for (auto left: lhs) {
        for (auto right: rhs) {
            result.push_back(left / right);
        }
    }

    return result;
}

optional<double> Div::optimize() {
    auto lhs_opt = lhs_->optimize();
    auto rhs_opt = rhs_->optimize();
    if (lhs_opt && rhs_opt) {
        return lhs_opt.value() / rhs_opt.value();
    } else if (lhs_opt) {
        lhs_ = std::make_unique<Number>(lhs_opt.value());
    } else if (rhs_opt) {
        rhs_ = std::make_unique<Number>(rhs_opt.value());
    }
    return nullopt;
}

std::string Div::print() const {
    return fmt::format("({} / {})", lhs_->print(), rhs_->print());
}

void Div::clear() {
    lhs_->clear();
    rhs_->clear();
}

NumericValues Pow::eval(const Frame& frame, const Match& match) const {
    auto lhs = lhs_->eval(frame, match);
    auto rhs = rhs_->eval(frame, match);

    auto result = NumericValues();
    result.reserve(lhs.size() * rhs.size());
    for (auto left: lhs) {
        for (auto right: rhs) {
            result.push_back(pow(left, right));
        }
    }

    return result;
}

optional<double> Pow::optimize() {
    auto lhs_opt = lhs_->optimize();
    auto rhs_opt = rhs_->optimize();
    if (lhs_opt && rhs_opt) {
        return pow(lhs_opt.value(), rhs_opt.value());
    } else if (lhs_opt) {
        lhs_ = std::make_unique<Number>(lhs_opt.value());
    } else if (rhs_opt) {
        rhs_ = std::make_unique<Number>(rhs_opt.value());
    }
    return nullopt;
}

std::string Pow::print() const {
    return fmt::format("{} ^({})", lhs_->print(), rhs_->print());
}

void Pow::clear() {
    lhs_->clear();
    rhs_->clear();
}

NumericValues Neg::eval(const Frame& frame, const Match& match) const {
    auto result = ast_->eval(frame, match);
    for (size_t i=0; i<result.size(); i++) {
        result[i] = -result[i];
    }
    return result;
}

optional<double> Neg::optimize() {
    auto optimized = ast_->optimize();
    if (optimized) {
        return - optimized.value();
    } else {
        return nullopt;
    }
}

std::string Neg::print() const {
    return fmt::format("(-{})", ast_->print());
}

void Neg::clear() {
    ast_->clear();
}

NumericValues Mod::eval(const Frame& frame, const Match& match) const {
    auto lhs = lhs_->eval(frame, match);
    auto rhs = rhs_->eval(frame, match);

    auto result = NumericValues();
    result.reserve(lhs.size() * rhs.size());
    for (auto left: lhs) {
        for (auto right: rhs) {
            result.push_back(fmod(left, right));
        }
    }

    return result;
}

optional<double> Mod::optimize() {
    auto lhs_opt = lhs_->optimize();
    auto rhs_opt = rhs_->optimize();
    if (lhs_opt && rhs_opt) {
        return fmod(lhs_opt.value(), rhs_opt.value());
    } else if (lhs_opt) {
        lhs_ = std::make_unique<Number>(lhs_opt.value());
    } else if (rhs_opt) {
        rhs_ = std::make_unique<Number>(rhs_opt.value());
    }
    return nullopt;
}

std::string Mod::print() const {
    return fmt::format("({} % {})", lhs_->print(), rhs_->print());
}

void Mod::clear() {
    lhs_->clear();
    rhs_->clear();
}

NumericValues Function::eval(const Frame& frame, const Match& match) const {
    auto result = ast_->eval(frame, match);
    for (size_t i=0; i<result.size(); i++) {
        result[i] = fn_(result[i]);
    }
    return result;
}

optional<double> Function::optimize() {
    auto optimized = ast_->optimize();
    if (optimized) {
        return fn_(optimized.value());
    } else {
        return nullopt;
    }
}

std::string Function::print() const {
    return fmt::format("{}({})", name_, ast_->print());
}

void Function::clear() {
    ast_->clear();
}

NumericValues Number::eval(const Frame& /*unused*/, const Match& /*unused*/) const {
    return NumericValues(value_);
}

optional<double> Number::optimize() {
    return value_;
}

std::string Number::print() const {
    if (std::round(value_) == value_) {
        return std::to_string(std::lround(value_));
    } else {
        return std::to_string(value_);
    }
}

NumericValues Distance::eval(const Frame& frame, const Match& match) const {
    auto results = NumericValues();
    for (auto i: i_.eval(frame, match)) {
        for (auto j: j_.eval(frame, match)) {
            results.push_back(frame.distance(i, j));
        }
    }
    return results;
}

std::string Distance::print() const {
    return fmt::format("distance({}, {})", i_.print(), j_.print());
}

NumericValues selections::Angle::eval(const Frame& frame, const Match& match) const {
    auto results = NumericValues();
    for (auto i: i_.eval(frame, match)) {
        for (auto j: j_.eval(frame, match)) {
            for (auto k: k_.eval(frame, match)) {
                results.push_back(frame.angle(i, j, k));
            }
        }
    }
    return results;
}

std::string selections::Angle::print() const {
    return fmt::format("angle({}, {}, {})", i_.print(), j_.print(), k_.print());
}

NumericValues selections::Dihedral::eval(const Frame& frame, const Match& match) const {
    auto results = NumericValues();
    for (auto i: i_.eval(frame, match)) {
        for (auto j: j_.eval(frame, match)) {
            for (auto k: k_.eval(frame, match)) {
                for (auto m: m_.eval(frame, match)) {
                    results.push_back(frame.dihedral(i, j, k, m));
                }
            }
        }
    }
    return results;
}

std::string selections::Dihedral::print() const {
    return fmt::format("dihedral({}, {}, {}, {})", i_.print(), j_.print(), k_.print(), m_.print());
}

NumericValues OutOfPlane::eval(const Frame& frame, const Match& match) const {
    auto results = NumericValues();
    for (auto i: i_.eval(frame, match)) {
        for (auto j: j_.eval(frame, match)) {
            for (auto k: k_.eval(frame, match)) {
                for (auto m: m_.eval(frame, match)) {
                    results.push_back(frame.out_of_plane(i, j, k, m));
                }
            }
        }
    }
    return results;
}

std::string OutOfPlane::print() const {
    return fmt::format("out_of_plane({}, {}, {}, {})", i_.print(), j_.print(), k_.print(), m_.print());
}

NumericValues NumericSelector::eval(const Frame& frame, const Match& match) const {
    return NumericValues(this->value(frame, match[argument_]));
}

optional<double> NumericSelector::optimize() {
    return nullopt;
}

std::string NumericSelector::print() const {
    return fmt::format("{}(#{})", name(), argument_ + 1);
}

std::string NumericProperty::name() const {
    if (is_ident(property_)) {
        return "[" + property_ + "]";
    } else {
        return "[\"" + property_ + "\"]";
    }
}

double NumericProperty::value(const Frame& frame, size_t i) const {
    const auto& property = frame[i].get(property_);
    if (property) {
        if (property->kind() == Property::DOUBLE) {
            return property->as_double();
        } else {
            throw selection_error(
                "invalid type for property [{}] on atom {}: expected double, got {}",
                property_, i, kind_as_string(property->kind())
            );
        }
    } else {
        auto residue = frame.topology().residue_for_atom(i);
        if (residue) {
            const auto& resProperty = residue->get(property_);
            if(!resProperty) {
                return nan("");
            }
            if (resProperty->kind() == Property::DOUBLE) {
                return resProperty->as_double();
            } else {
                throw selection_error(
                    "invalid type for property [{}] on the residue containing atom {}: expected double, got {}",
                    property_, i, kind_as_string(resProperty->kind())
                );
            }
        }
        // No property with the given name
        // return nan so that all comparaison down the line evaluate to false
        return nan("");
    }
}

std::string Index::name() const {
    return "index";
}

double Index::value(const Frame& /*unused*/, size_t i) const {
    return static_cast<double>(i);
}

std::string Resid::name() const {
    return "resid";
}

double Resid::value(const Frame& frame, size_t i) const {
    auto residue = frame.topology().residue_for_atom(i);
    if (residue && residue->id()) {
        return static_cast<double>(*residue->id());
    } else {
        return -1;
    }
}

std::string Mass::name() const {
    return "mass";
}

double Mass::value(const Frame& frame, size_t i) const {
    return frame[i].mass();
}

std::string Position::name() const {
    switch (coordinate_) {
    case Coordinate::X:
        return "x";
    case Coordinate::Y:
        return "y";
    case Coordinate::Z:
        return "z";
    default:
        unreachable();
    }
}

double Position::value(const Frame& frame, size_t i) const {
    return frame.positions()[i][static_cast<size_t>(coordinate_)];
}

std::string Velocity::name() const {
    switch (coordinate_) {
    case Coordinate::X:
        return "vx";
    case Coordinate::Y:
        return "vy";
    case Coordinate::Z:
        return "vz";
    default:
        unreachable();
    }
}

double Velocity::value(const Frame& frame, size_t i) const {
    if (frame.velocities()) {
        const auto& velocities = *frame.velocities();
        return velocities[i][static_cast<size_t>(coordinate_)];
    } else {
        // return nan so that all comparaison down the line evaluate to false
        return std::nan("");
    }
}

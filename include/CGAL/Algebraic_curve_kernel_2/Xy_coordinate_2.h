// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     : Eric Berberich <eric@mpi-inf.mpg.de>
//                 Pavel Emeliyanenko <asm@mpi-sb.mpg.de>
//
// ============================================================================

#ifndef CGAL_ALGEBRAIC_CURVE_KERNEL_XY_COORDINATE_2_H
#define CGAL_ALGEBRAIC_CURVE_KERNEL_XY_COORDINATE_2_H

#include <CGAL/basic.h>
#include <boost/pool/pool_alloc.hpp>
#include <boost/numeric/interval.hpp>

#include <CGAL/Bbox_2.h>

CGAL_BEGIN_NAMESPACE

namespace CGALi {

template < class AlgebraicCurveKernel_2, class Rep_, 
      class HandlePolicy_ ,
      class Allocator_>
        //::boost::fast_pool_allocator<Rep_> >
class Xy_coordinate_2;


template < class AlgebraicCurveKernel_2 >
class Xy_coordinate_2_rep {

public:
    // this first template argument
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;

    // myself
    typedef Xy_coordinate_2_rep<Algebraic_curve_kernel_2> Self;

    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
        Curve_analysis_2;

    typedef typename Algebraic_curve_kernel_2::Algebraic_real_1 
        Algebraic_real_1;

    typedef CGAL::Bbox_2 Bbox_2;

    typedef CGAL::Handle_with_policy<Self>
        Xy_coordinate_2_inst;

    // constructors
public:
    // default constructor ()
    Xy_coordinate_2_rep() : _m_arcno(-1) {
    }
    
    // standard constructor
    Xy_coordinate_2_rep(const Algebraic_real_1& x,
        const Curve_analysis_2& curve, int arcno) :
            _m_x(x), _m_curve(curve), _m_arcno(arcno) {
    }

    // data
    // x-coordinate
    Algebraic_real_1 _m_x;
    
    // supporting curve
    mutable Curve_analysis_2 _m_curve;
    
    // arc number on curve
    mutable int _m_arcno;

    // y-coordinate
    mutable boost::optional< Algebraic_real_1 > _m_y;

    //! A bounding box for the given point
    mutable boost::optional< std::pair<double,Bbox_2> > _m_bbox_2_pair;

};

//! \brief class \c Xy_coordinate_2 represents a single root of a system of 
//! two polynomial equations in two variables that are models 
//! \c AlgebraicCurveKernel_2::Polynomial_2
//!
//! \c Xy_coordinate_2 coordinate is represented by an \c Algebraic_real_1,
//! a supporting curve and an arcno and is valid only for finite solutions,
//! i.e., it cannot represent points at infinity 
template <class AlgebraicCurveKernel_2, 
          class Rep_ = CGALi::Xy_coordinate_2_rep<AlgebraicCurveKernel_2>,
          class HandlePolicy_= CGAL::Handle_policy_union, 
          class Allocator_ = CGAL_ALLOCATOR(Rep_) >
class Xy_coordinate_2 : 
    public ::CGAL::Handle_with_policy<Rep_, HandlePolicy_, Allocator_> 
{
public:
    //! \name public typedefs
    //!@{
    
    //! this instance's first template parameter
    typedef AlgebraicCurveKernel_2 Algebraic_curve_kernel_2;
    
    //! this instance's second template parameter
    typedef Rep_ Rep;
    
    //! this instance's third template parameter
    typedef HandlePolicy_ Handle_policy;
    
    //! this instance's fourth template parameter
    typedef Allocator_ Allocator;

    //! this instance itself
    typedef Xy_coordinate_2<Algebraic_curve_kernel_2, Rep, Handle_policy,
        Allocator> Self;
        
    //! an instance of AlgebraicKernel_1
    typedef typename Algebraic_curve_kernel_2::Algebraic_kernel_1 
        Algebraic_kernel_1;
    
    //! type of (explicit) x- and y-coordinates
    typedef typename Algebraic_curve_kernel_2::Algebraic_real_1 
        Algebraic_real_1;

    //! Coefficient type
    typedef typename Algebraic_curve_kernel_2::Coefficient Coefficient;

    //! type of curve pair analysis
    typedef typename Algebraic_curve_kernel_2::Curve_pair_analysis_2
                Curve_pair_analysis_2;
    
    //! type of curve analysis
    typedef typename Algebraic_curve_kernel_2::Curve_analysis_2
                Curve_analysis_2;
    
    //! the handle superclass
    typedef ::CGAL::Handle_with_policy<Rep, Handle_policy, Allocator> Base;

    typedef typename Algebraic_curve_kernel_2::X_real_traits_1 X_real_traits_1;

    //! Algebraic_real_traits for this type
    typedef typename Algebraic_curve_kernel_2::Y_real_traits_1 Y_real_traits_1;

    //! type for approximation boundaries
    typedef typename X_real_traits_1::Bound Bound;

    //! type for bound intervals
    typedef boost::numeric::interval<Bound> Bound_interval;

    typedef CGAL::Coercion_traits<Coefficient,Bound> Coercion;

    typedef typename CGAL::Coercion_traits<Coefficient,Bound>::Type
        Coercion_type;

    typedef boost::numeric::interval<Coercion_type> Coercion_interval;

    //! Type for the bounding box
    typedef typename Rep::Bbox_2 Bbox_2;
    
    //!@}
private:
    //! \name private methods
    //!@{

    /*!\brief
     * Simplifies the representation of two points whose supporting curves
     * share a common part.
     */
    static bool _simplify(const Xy_coordinate_2& p, const Xy_coordinate_2& q) 
    {
        std::vector<Curve_analysis_2> parts_of_f, parts_of_g, common;
        Algebraic_curve_kernel_2 ak_2;

        if(ak_2.decompose_2_object()(p.curve(), q.curve(), 
            std::back_inserter(parts_of_f), std::back_inserter(parts_of_g),
                std::back_inserter(common))) {

            CGAL_assertion((parts_of_f.size() == 1 ||
                       parts_of_g.size() == 1) && common.size() == 1);
            if(parts_of_f.size() == 1) {
                p.simplify_by(ak_2.construct_curve_pair_2_object()(
                    parts_of_f[0], common[0]));
            } 
            
            if(parts_of_g.size() == 1) {
                q.simplify_by(ak_2.construct_curve_pair_2_object()(
                    parts_of_g[0], common[0]));
            } 
            return true;
        }
        return false;
    }
    
    //!@}
public:
    //!\name Constructors
    //!@{

    /*!\brief 
     * default constructor
     *
     * A default-constructed point supports no operation other than
     * having \c CGAL::degree(curve()) return \c -1. 
     */
    Xy_coordinate_2() : 
        Base(Rep()) { 
    }

    /*!\brief
     * copy constructor
     */
    Xy_coordinate_2(const Self& p) : 
        Base(static_cast<const Base&>(p)) {  
    }

    /*!\brief
     * Point at \c x, on \c curve with \c arcno. Finite points on vertical arcs
     * are also constructed in this way
     */
     // TODO: construct this from curve analysis object ?
    Xy_coordinate_2(const Algebraic_real_1& x, const Curve_analysis_2& curve,
                 int arcno) :
            Base(Rep(x, curve, arcno)) {
            
        CGAL_precondition(arcno >= 0);
        CGAL_precondition_code(
            typename Curve_analysis_2::Status_line_1 v =
                curve.status_line_for_x(x);
        );
        CGAL_precondition(arcno >= 0 && arcno < v.number_of_events());
    }
    
    /*!\brief
     * constructs a point from a given represenation
     */
    Xy_coordinate_2(Rep rep) : 
        Base(rep) {  
    }
   
    //!@}
public:
    //!\name Access functions
    //!@{
    
    /*!\brief 
     * x-coordinate of the point
     */
    const Algebraic_real_1& x() const { 
        return this->ptr()->_m_x; 
    }

    /*!
     * \brief y-coordinate of this point
     *
     * Note: In general, this method results in a extremly large polynomial
     * for the y-coordinate. It is recommended to use it carefully,
     * and using get_approximation_y() instead whenever approximations suffice.
     */
    Algebraic_real_1 y() const {

        typedef typename Algebraic_curve_kernel_2::Polynomial_2 Polynomial_2;
        typedef typename Algebraic_curve_kernel_2::Polynomial_1 Polynomial_1;
        
        typedef std::vector< Algebraic_real_1 > Roots;
        typedef typename Curve_analysis_2::Status_line_1 Key;
        typedef Roots Data;
        typedef std::map< Key, Data, CGAL::Handle_id_less_than< Key > > 
            Y_root_map;
        
        static Y_root_map y_root_map;

        if (!this->ptr()->_m_y) {
            
            Polynomial_2 f = curve().primitive_polynomial_2();
            // This will be the defining polynomial of y
            Polynomial_1 y_pol;

            // Filter: If we know that the point is critical, we can use
            // the resultant of f and f_y with respect to x as polynomial
            bool point_is_certainly_critical = false;
            typename Curve_analysis_2::Status_line_1 line =
                curve().status_line_at_exact_x(x());
            
            typename Y_root_map::iterator yit = 
                y_root_map.find(line);

            // TODO: Cache resultant computation
            // exacus-related code shouldn't be used here
            //curve().x_to_index(x(),i,is_event);
            if (line.is_event()) {
                //typename Internal_curve_2::Event1_info ev_info =
                //   curve().event_info(i);
                typename Curve_analysis_2::Status_line_1::Arc_pair ipair =
                    line.number_of_incident_branches(arcno());
                
                if (ipair.first != 1 || ipair.second != 1) {
                    point_is_certainly_critical = true;
                    y_pol = 
                        typename CGAL::Polynomial_traits_d< Polynomial_1 >::Make_square_free()(
                                CGAL::CGALi::resultant
                                (typename CGAL::Polynomial_traits_d
                                 <Polynomial_2>::Swap() (f,0,1),
                                 typename CGAL::Polynomial_traits_d
                                 <Polynomial_2>::Swap() 
                                 (CGAL::diff(f),0,1))
                        );
                    // BUGFIX: y_pol might be zero:
                    if(y_pol.is_zero()) {
                        // force re-computation with bigger resultant
                        point_is_certainly_critical=false;
                    }                             
                    
                }
            }
            
            if (!point_is_certainly_critical) {
                
                Polynomial_2 r(x().polynomial());
                y_pol = typename CGAL::Polynomial_traits_d< Polynomial_1 >::Make_square_free()(
                        CGAL::CGALi::resultant
                        (typename CGAL::Polynomial_traits_d
                         <Polynomial_2>::Swap() (f,0,1),
                         typename CGAL::Polynomial_traits_d
                         <Polynomial_2>::Swap() (r,0,1))
                );
                
            }
            typename Algebraic_kernel_1::Solve_1 real_roots;
            
            Roots y_roots;
            real_roots(y_pol, std::back_inserter(y_roots), false ); 
            
            Bound_interval y_iv = get_approximation_y();
            
            typedef typename std::vector<Algebraic_real_1>::const_iterator
                Iterator;
            
            std::list< Iterator > candidates;
            
            for (Iterator it = y_roots.begin(); it != y_roots.end(); it++) {
                Bound_interval it_interval(it->low(), it->high());
                if (boost::numeric::overlap(it_interval, y_iv)) {
                    candidates.push_back(it);
                }
            }
            CGAL_assertion(!candidates.empty());

            while (candidates.size() > 1) {
                refine_y();
                y_iv = get_approximation_y();

                for (typename std::list< Iterator >::iterator dit, cit =
                         candidates.begin(); cit != candidates.end(); ) {
                    bool remove = false;
                    Bound_interval 
                        cit_interval((*cit)->low(), (*cit)->high());
                    if (!boost::numeric::overlap(cit_interval, y_iv)) {
                        dit = cit;
                        remove = true;
                    }
                    cit++;
                    if (remove) {
                        candidates.erase(dit);
                    }
                }
            }
            CGAL_assertion(static_cast< int >(candidates.size()) == 1);
            this->ptr()->_m_y = 
                Algebraic_real_1(
                        (*candidates.begin())->polynomial(), 
                        (*candidates.begin())->low(), 
                        (*candidates.begin())->high()
                );
        }
        CGAL_postcondition(this->ptr()->_m_y);
        return *this->ptr()->_m_y;
    }
    
    /*!\brief
     * supporting curve of the point
     */
    Curve_analysis_2 curve() const {
        return this->ptr()->_m_curve; 
    }
    
    /*!\brief
     * arc number of point
     *
     */
    int arcno() const { 
        return this->ptr()->_m_arcno; 
    }

    //!@}
public:
    //!\name comparison predicates
    //!@{

    /*!\brief
     * compares x-coordinates of \c *this with \c q
     * 
     * do we need this method or one should use Algebraic_curve_kernel_2
     * directly ?
     */
    CGAL::Comparison_result compare_x(const Self& q) const {
    
        if(this->is_identical(q)) {
            return CGAL::EQUAL;
        }
        Algebraic_curve_kernel_2 ak_2;
        return ak_2.compare_x_2_object()(this->x(), q.x());
    }

    /*!\brief
     * compares \c *this with \c q lexicographically
     */
    CGAL::Comparison_result compare_xy(const Self& q, 
        bool equal_x = false) const {
        
        if(is_identical(q)) 
            return CGAL::EQUAL;

        CGAL::Comparison_result res = (equal_x ? CGAL::EQUAL : compare_x(q)); 
        if(res == CGAL::EQUAL) {
            res = _compare_y_at_x(q);
        }
        return res;
    }
    
    //! equality
    bool operator == (const Self& q) const {return q.compare_xy(*this)== 0;}
    
    //! inequality
    bool operator != (const Self& q) const {return q.compare_xy(*this)!= 0;}

    //! less than in (x,y) lexicographic order
    bool operator <  (const Self& q) const {return q.compare_xy(*this)> 0;}

    //! less-equal in (x,y) lexicographic order
    bool operator <= (const Self& q) const {return q.compare_xy(*this)>= 0;}

    //! greater than in (x,y) lexicographic order
    bool operator >  (const Self& q) const {return q.compare_xy(*this)< 0;}

    //! greater-equal in (x,y) lexicographic order
    bool operator >= (const Self& q) const {return q.compare_xy(*this)<= 0;}
    
    //!@}
private:
    //!@{
    //! \name 

    /*!\brief
     * compares y-coordinates for covertical points \c *this and \c q
     *
     * \pre x() == q.x()
     */
    CGAL::Comparison_result _compare_y_at_x(const Self& q) const 
    {
        CGAL_precondition(this->compare_x(q) == CGAL::EQUAL);
    
        Curve_analysis_2 f = curve(), g = q.curve();
        if(f.is_identical(g)) 
            return CGAL::sign(arcno() - q.arcno());
        if(Self::_simplify(*this, q)) 
            // restart since supporting curves might be equal now
            return _compare_y_at_x(q);
                        
        Algebraic_curve_kernel_2 ak_2;
        Curve_pair_analysis_2 cpa_2 =
            ak_2.construct_curve_pair_2_object()(f, g);
            
            
        typename Curve_pair_analysis_2::Status_line_1 vline =
            cpa_2.status_line_for_x(x());
        return CGAL::sign(vline.event_of_curve(arcno(), f) -
                    vline.event_of_curve(q.arcno(), g));
    }
    
    //!@}
public:
    //!\name Reconstructing functions
    //!@{
    
    /*!\brief
     * Simplifies the representation of a point.
     * 
     * Given a decomposition of the point's supporting \c curve() into 
     * a pair of two curves \c pair, this function searches this point
     * in the curve pair and resets the curve and the arcno to this
     * found arc. It can happen, that both curves of the pair fit this 
     * condition (intersection of the two curves at this point), then it
     * chooses the simpler one (less total degree).
     *
     * \pre pair must be a decomposition of curve()
     */
    void simplify_by(const Curve_pair_analysis_2& cpa_2) const {
    
        typedef typename Algebraic_curve_kernel_2::Polynomial_2 Poly_2;
        
        CGAL_precondition_code(
            typename Curve_analysis_2::Polynomial_2 mult =
                    cpa_2.curve_analysis(0).polynomial_2() *
                    cpa_2.curve_analysis(1).polynomial_2();
        );
        // common parts
        CGAL_precondition(CGAL::CGALi::resultant(mult,
                                                 curve().polynomial_2())
                          .is_zero());
        // full parts
        CGAL_precondition(CGAL::degree(mult) == 
                          CGAL::degree(curve().polynomial_2()));
        CGAL_precondition(CGAL::total_degree(mult) ==
                          CGAL::total_degree(curve().polynomial_2()));

        typename Curve_pair_analysis_2::Status_line_1 cpv_line =
                cpa_2.status_line_for_x(x());
        // # of arcs must match
        CGAL_precondition_code(
            typename Curve_analysis_2::Status_line_1 cv_line =
                curve().status_line_for_x(x());
        );
        CGAL_precondition(cpv_line.number_of_events() == 
            cv_line.number_of_events());

        int cid = 0;
        std::pair<int, int> p = cpv_line.curves_at_event(arcno());
        if(p.first != -1 && p.second != -1) {
            // both curves involved: choose simpler one
            // Remark: In this case, a vertical line in the curves can be
            // ignored, since it has not been considered when constructing
            // the point from the composed curved (also including this vertical
            // line). Therefore, the old arc number is also valid in the curve
            // pair.
            Poly_2 ff = cpa_2.curve_analysis(0).polynomial_2(),
                   gg = cpa_2.curve_analysis(1).polynomial_2();
            if(total_degree(ff) > total_degree(gg)) 
                cid = 1;
        } else 
            cid = (p.first != -1 ? 0 : 1);
        // overwrite data
        this->ptr()->_m_curve = cpa_2.curve_analysis(cid);
        this->ptr()->_m_arcno = (cid == 0 ? p.first : p.second);
    }
    
    //! befriending output iterator
   // friend std::ostream& operator << <>(std::ostream& os, const Self& pt);

    //!@}
public:
    
    //! Returns whether the x-coordinate equals zero
    bool is_x_zero() const {
        return this->ptr()->_m_x.is_zero();
    }

    //! Returns whether the y-coordinate equals zero
    bool is_y_zero() const {
        Bound_interval y_iv = get_approximation_y();
        if( y_iv.lower() > 0 || y_iv.upper() < 0 ) {
            return false;
        }
        CGAL_assertion(CGAL::degree(curve().polynomial_2())>=0);
        
        typename Curve_analysis_2::Polynomial_2::NT constant_pol
            = curve().polynomial_2()[0];
        bool zero_is_root_of_local_pol 
            = this->ptr()->_m_x.is_root_of(constant_pol);
        // Since we know that y_iv is an _isolating_ interval,
        // we can immediately return
        return zero_is_root_of_local_pol;
        
    }
    
    // returns a double approximation of the point
    std::pair<double, double> to_double() const {

        typedef typename Get_arithmetic_kernel<Bound>::Arithmetic_kernel AT;
        typedef typename AT::Bigfloat_interval BFI; 
        typedef typename CGAL::Bigfloat_interval_traits<BFI>::Bound BF;

        long old_prec = get_precision(BFI());
        
        set_precision (BFI(), 53);

        double double_x = this->ptr()->_m_x.to_double();
        double double_y;

        typename Y_real_traits_1::Lower_bound lower;
        typename Y_real_traits_1::Upper_bound upper;
        typename Y_real_traits_1::Refine refine;

        if (lower(*this)==upper(*this)) {
            double_y = CGAL::to_double(convert_to_bfi(lower(*this)));
        } else if(is_y_zero()) {
            double_y = 0.;
        } else {
            while(CGAL::sign(lower(*this)) != 
                  CGAL::sign(upper(*this)) ) {
                refine(*this);
            }
            long final_prec = set_precision(BFI(),get_precision(BFI())+4);
            
            BFI bfi = CGAL::hull(convert_to_bfi(lower(*this)), 
                                  convert_to_bfi(upper(*this)));
            
            while( !singleton(bfi) &&  
                   get_significant_bits(bfi) < final_prec  ){
                refine(*this);
                bfi = CGAL::hull(
                        convert_to_bfi(lower(*this)), 
                        convert_to_bfi(upper(*this)));
            }
            double_y 
                = CGAL::to_double((CGAL::lower(bfi)+ CGAL::upper(bfi)) / 2);
        }
        set_precision(BFI(),old_prec);
        return std::make_pair(double_x, double_y); 
    }

    public:
    //!\name Approximating functions
    //!@{

    /*!
     * \brief gets approximation of x
     */
    Bound_interval get_approximation_x() const {
        
        typename X_real_traits_1::Lower_bound lower;
        typename X_real_traits_1::Upper_bound upper;

        return Bound_interval(lower(this->ptr()->_m_x), 
                                 upper(this->ptr()->_m_x));
        
    }

    /*!
     * \brief gets approximation of x that is smaller than bound
     */
    Bound_interval get_approximation_x(Bound bound) const {
        
        CGAL_assertion(bound > 0);

        typename X_real_traits_1::Lower_bound lower;
        typename X_real_traits_1::Upper_bound upper;
        typename X_real_traits_1::Refine refine;

        while(upper(this->ptr()->_m_x) - lower(this->ptr()->_m_x) >= bound) {
            refine(this->ptr()->_m_x);
        }
        return Bound_interval(lower(this->ptr()->_m_x), 
                                 upper(this->ptr()->_m_x));
    
    }

    /*!
     * \brief gets approximation of y
     *
     */
    Bound_interval get_approximation_y() const {
        typename Y_real_traits_1::Lower_bound lower;
        typename Y_real_traits_1::Upper_bound upper;
        return Bound_interval(lower(*this), upper(*this));
    }

    /*!
     * \brief gets approximation of y that is smaller than bound
     */
    Bound_interval get_approximation_y(Bound bound) const {
        
        CGAL_assertion(bound > 0);

        typename Y_real_traits_1::Lower_bound lower;
        typename Y_real_traits_1::Upper_bound upper;
        typename Y_real_traits_1::Refine refine;

        while(upper(*this) - lower(*this) >= bound) {
            refine(*this);
        }
        return Bound_interval(lower(*this),upper(*this));    
    }
    
    /*!\brief 
     *  computes at least half's the current interval of the first coordinate 
     */
    void refine_x() const {
        // typename Algebraic_kernel_1::Refine refine
        // refine(this->ptr()->_m_x)
        //this->ptr()->_m_x.refine();
        typename X_real_traits_1::Refine()(this->ptr()->_m_x);
    }
    
    /*!\brief 
     *  refines the first coordinate w.r.t relative precision
     */
    void refine_x(int rel_prec) {
        // typename Algebraic_kernel_1::Refine refine
        // refine(this->ptr()->_m_x)
        typename X_real_traits_1::Refine()(this->ptr()->_m_x, rel_prec);
    }

    /*!
     * \brief Refines the x-xoordinate
     */
    void refine_y() const {
        typename Y_real_traits_1::Refine refine;
        refine(*this);
    }

    Bbox_2 approximation_box_2() const {
           
        double x_min, x_max, y_min, y_max;
        
        typename X_real_traits_1::Lower_bound x_lower;
        typename X_real_traits_1::Upper_bound x_upper;
        
        typename Y_real_traits_1::Lower_bound y_lower;
        typename Y_real_traits_1::Upper_bound y_upper;
        
        x_min = CGAL::to_interval(x_lower(this->x())).first;
        x_max = CGAL::to_interval(x_upper(this->x())).second;
        y_min = CGAL::to_interval(y_lower(*this)).first;
        y_max = CGAL::to_interval(y_upper(*this)).second;
        return Bbox_2(x_min,y_min,x_max,y_max);
    }

    Bbox_2 approximation_box_2(double b) const {

        typename X_real_traits_1::Refine x_refine;
        typename Y_real_traits_1::Refine y_refine;

        if(this->ptr()->_m_bbox_2_pair) {
            double cached_prec = this->ptr()->_m_bbox_2_pair.get().first;
            if(cached_prec <= b) {
                return this->ptr()->_m_bbox_2_pair.get().second;
            }
        }

        Bbox_2 box = approximation_box_2();
        double x_diff = box.xmax()-box.xmin();
        double y_diff = box.ymax()-box.ymin();
        while(x_diff >= b || y_diff >= b) {
            if(x_diff >= y_diff) {
                x_refine(this->x());
            } else {
                y_refine(*this);
            }
            box = approximation_box_2();
            x_diff = box.xmax()-box.xmin();
            y_diff = box.ymax()-box.ymin();
        }
        this->ptr()->_m_bbox_2_pair = std::make_pair(b,box);
        return box;
    }

    template<typename Polynomial_2>
    Coercion_interval interval_evaluate_2(const Polynomial_2& p) const {
        
        typename Coercion::Cast cast;

        typedef typename CGAL::Polynomial_traits_d<Polynomial_2>::
            Coefficient_const_iterator
            Coefficient_const_iterator;
        
        typedef typename CGAL::Polynomial_traits_d<Polynomial_2>
            ::Coefficient_const_iterator_range 
            Coefficient_const_iterator_range;
        
        Coercion_interval iy(cast(get_approximation_y().lower()),
                             cast(get_approximation_y().upper()));

        // CGAL::Polynomial does not provide Coercion_traits for number
        // types => therefore evaluate manually
        Coefficient_const_iterator_range range = 
            typename CGAL::Polynomial_traits_d<Polynomial_2>
            :: Construct_coefficient_const_iterator_range()(p);
        
        Coefficient_const_iterator it = range.second - 1;
        
        Coercion_interval res(interval_evaluate_1(*it));
        
        Coefficient_const_iterator p_begin = range.first;

        while((it--) != p_begin) 
            res = res * iy + (interval_evaluate_1(*it));
        return res;
    }
    
    template<typename Polynomial_1>
    Coercion_interval interval_evaluate_1(const Polynomial_1& p) const {
        
        typename Coercion::Cast cast;

        typedef typename CGAL::Polynomial_traits_d<Polynomial_1>
            ::Coefficient_const_iterator Coefficient_const_iterator;

        Coercion_interval ix(cast(get_approximation_x().lower()),
                             cast(get_approximation_x().upper()));
        
        typedef typename CGAL::Polynomial_traits_d<Polynomial_1>
            ::Coefficient_const_iterator_range 
            Coefficient_const_iterator_range;
        
        Coefficient_const_iterator_range range = 
            typename CGAL::Polynomial_traits_d<Polynomial_1>
            :: Construct_coefficient_const_iterator_range()(p);
        
        Coefficient_const_iterator it = range.second - 1;
        
        Coercion_interval res(cast(*it));

        Coefficient_const_iterator p_begin = range.first;
        while((it--) != p_begin) 
            res = res * ix + Coercion_interval(cast(*it));
        return res;
    }

     // friend function to provide a fast hashing
    friend std::size_t hash_value(const Self& x) {
        return static_cast<std::size_t>(x.id());
    }

    //!@}

}; // class Xy_coordinate_2

template < class AlgebraicCurveKernel_2, class Rep> 
std::ostream& operator<< (std::ostream& os, 
    const Xy_coordinate_2<AlgebraicCurveKernel_2, Rep>& pt)
{
    if(::CGAL::get_mode(os) == ::CGAL::IO::PRETTY) {
        os << "[x-coord: " << CGAL::to_double(pt.x()) << "; curve: " <<
            pt.curve().polynomial_2() << 
            "; arcno: " << pt.arcno() << "]\n";
    } else { // ASCII output
        os << "[x-coord: " << CGAL::to_double(pt.x()) << "; curve: " <<
            pt.curve().polynomial_2() << 
            "; arcno: " << pt.arcno() << "]\n";
    }
    return os;    
}

} // namespace CGALi

CGAL_END_NAMESPACE

#endif // CGAL_ALGEBRAIC_CURVE_KERNEL_XY_COORDINATE_2_H

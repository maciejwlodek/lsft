"""
Given the front projection of a Legendrian knot as the plat closure
of a braid, compute its classical invariants, Chekanov-Eliashberg
algebra, ruling polynomials, and augmentation category.

This code was originally written for the paper "A bordered
Chekanov-Eliashberg algebra", arXiv:1004.4929, and can also compute
the augmentation categories described in "Augmentations are sheaves"
(with Ng, Rutherford, Shende, and Zaslow), arXiv:1502.04939.

AUTHORS:

- Steven Sivek (2010-10-10, updated 2013-09-21 and 2014-11-05;
  major update for augmentation categories 2015-05-03)

The knot data was tabulated by Sumana Shrestha and taken from the
Legendrian Invariants notebook on Josh Sabloff's website, with the
exception of the m(9_36) plat.  Its introductory text reads:

"The following are Legendrian plat representations of maximal Thurston-
Bennequin Legendrian knots and their mirrors up to 9 crossings. There
is only one Legendrian representative per knot, even though a given 
knot may have many different Legendrian representatives."

- minor update by Maciej Wlodek for Python3 compatibility and to add plat plotting functionality


EXAMPLES::

   sage: K = plat([2,2,2])
   sage: K.gradings()
   (0, 0, 0, 1, 1)
   sage: K.differentials()
   (0, 0, 0, t + a0 + a2 + a0*a1*a2, 1 - a0 - a2 - a2*a1*a0)
   
   sage: K = LegendrianKnotTable['4_1']; K
   4_1
   sage: K.crossings()
   (2, 1, 4, 5, 1, 1, 3, 5, 1, 2, 5, 4)
   sage: K.tb(), K.rot()
   (-3, 0)
   
   sage: LegendrianKnotTable['m8_21'].all_homology()
   [1, t + 2]
   sage: LegendrianKnotTable['m8_21'].ruling_polynomial()
   3 + 2*z^2
   sage: LegendrianKnotTable['m8_21'].ncable(2).ruling_polynomial()
   9*z^-1 + 41*z + 70*z^3 + 47*z^5 + 12*z^7 + z^9

   sage: K = LegendrianKnotTable['m9_42']; K.all_homology()
   []
   sage: W = K.whitehead_double()
   sage: a = map(GF(2), [1,1] + [1 if i in [7,14,19,23,30,34,50,52,57,60] else 0 for i in range(69)])
   sage: W.homology(a)
   t^2 + t + 1

   sage: A = LegendrianKnotTable['m3_1'].augcat(GF(2)); A
   Aug(m3_1, GF(2))
   sage: A.knot().ruling_polynomial()
   2 + z^2
   sage: A.cardinality()
   5
"""

class plat:
    """
    Plat diagrams.  These are specified by a list of positive integers
    representing a positive braid whose plat closure is the desired front.
    All fronts are assumed to be knot diagrams, i.e. they have only one
    component.

    Note that plat closure is not the more common notion of braid closure:
    we assume there are an even number 2n of strands, and then we add n left
    cusps and n right cusps to the braid joining strands (2i-1, 2i) on
    either end for 1 <= i <= n.  For example, the plat closure of [2,2,2]
    has 4 strands, with three crossings between strands 2 and 3, and two
    left cusps and two right cusps; it is a Legendrian trefoil.

    EXAMPLES::

       sage: K = LegendrianKnotTable['m3_1']; K.crossings()
       (2, 2, 2)
       sage: K.tb(), K.rot()
       (1, 0)
       sage: len(K.augmentations())
       5
       sage: K.all_homology()
       [1]
       sage: K.augmentation_category(GF(3)).cardinality()
       5
       sage: W = K.whitehead_double()
       sage: W.all_homology(reduced = False)
       [t + 2, 3*t + 6 + 2*t^-1, 2*t + 4 + t^-1]
    """

    def __init__(self, crossinglist, name=None, **extra_args):
        """
        Create a plat with n left cusps.  The crossings are specified
        from left to right by integers between 1 and 2n-1, inclusive,
        where i denotes the ith strand crossing over the (i+1)st strand.
        Strands are numbered from 1 at the top to 2n at the bottom.

        INPUT:

        - ``crossinglist`` -- a list of positive integers describing the
          crossings from left to right.

        - ``name`` -- an optional string such as '3_1'.

        EXAMPLES::

           sage: K = plat([2,2,2]) # a Legendrian trefoil (the mirror of 3_1)
        """
        self.__crossinglist__ = tuple(crossinglist)
        if len(crossinglist) == 0:
            self.__ncusps__ = 1
        else:
            self.__ncusps__ = ceil((max(crossinglist)+1)/2)
        n = len(self.__crossinglist__) + self.__ncusps__
        varnames = ['t','u'] 
        varnames += ['a0'] if n==1 else ['a'+str(i) for i in range(n)]
        self.__ring__ = FreeAlgebra(ZZ, n+2,varnames)
        self.__ab_ring__ = PolynomialRing(ZZ, n+2,varnames)
        self.__augmentations__ = {}
        self.__all_homology__ = {}
        self.__ncopy_differential_list__ = {}
        self.__name__ = copy(name)
        self.__rp__ = {}
        self.__augcats__ = {}
        self.__ncable__ = 1
        
        # Orientation of top strand at the top right cusp
        self.__top_oriented_left_to_right__ = True
        if 'top_oriented_l2r' in extra_args.keys():
            self.__top_oriented_left_to_right__ = extra_args['top_oriented_l2r']

        # Compute gradings -- this lets us check that we have an actual
        # knot instead of a link with several components
        _ = self.gradings()

    def __repr__(self):
        """
        Return a string describing this plat.

        EXAMPLES::

           sage: plat([2,2,2,2,2])
           plat([2, 2, 2, 2, 2])
           sage: plat([2,2,2,2,2], 'cinquefoil')
           cinquefoil
        """
        if self.name() is None:
          return "plat(%s)"%(list(self.crossings()))
        return self.name()

    def __abelian_differentials__(self):
        """
        Return the list of differentials in the abelianization of the DGA.
        """
        try:
            return self.__comm_d_list__
        except AttributeError:
            R = self.__comm_ring__()
            self.__comm_d_list__ = tuple([R(d) for d in self.differentials()])
            return self.__comm_d_list__

    def __comm_ring__(self):
        """
        Return the commutative ring over which the abelianized DGA is defined.
        """
        return self.__ab_ring__

    def all_homology(self, base_ring = GF(2), rho = 0,
                     reduced = True, verbose = False):
        """
        Return the list of all Poincare polynomials of rho-graded
        augmentations of K.  If reduced is True, return the reduced
        polynomials p(t) (i.e. the actual polynomials for which
        t+p(t)+p(1/t) is the Poincare polynomial), and return the
        Poincare polynomial otherwise.

        INPUT:

        - ``base_ring`` -- (default: GF(2)) A finite field.

        - ``rho`` -- (default: 0) A nonnegative integer.
        
        - ``reduced`` -- (default: True) a boolean flag indicating that we
          should return polynomials p(t) corresponding to Poincare
          polynomials t+p(t)+p^(t^(-1)); if False, we return the entire
          Poincare polynomial.

        - ``verbose`` -- (default: False) boolean

        OUTPUT:

        A list of polynomials in ZZ[t], or Laurent polynomials if reduced
        is False.

        EXAMPLES::

           sage: K = plat([4,3,2,1,1,5,3,2,5,4,4,4,3,2,5,4])
           sage: K.all_homology()
           [1, t + 2]
           sage: K.all_homology(reduced = False)
           [t + 2, 2*t + 4 + t^-1]
        """
        try:
            hset = self.__all_homology__[(base_ring, rho)]
        except KeyError:
            if self.rot() != 0:
                hset = Set([])
            else:
                elist = self.augmentations(base_ring, rho, verbose)
                homlist = []
                for e in elist:
                    h = self.homology(e, False)
                    if h not in homlist:
                        homlist.append(h)
                hset = Set(homlist)
            self.__all_homology__[(base_ring, rho)] = hset

        if reduced:
            fn = lambda p: p
        else:
            S.<t> = LaurentPolynomialRing(ZZ, 't')
            fn = lambda p: t + p(t) + p(t^(-1))
        return list(map(fn, hset))

    def augcat(self, base_ring = GF(2), rho = 0):
        """
        A synonym for augmentation_category.
        """
        return self.augmentation_category(base_ring, rho)
    
    def augmentation_category(self, base_ring = GF(2), rho = 0):
        """
        Return the rho-graded augmentation category of this plat.

        INPUT:

        - ``base_ring`` -- (default: GF(2)) A field

        - ``rho`` -- (default: 0) A nonnegative integer

        OUTPUT:

        An instance of class AugmentationCategory

        EXAMPLES::

           sage: A = LegendrianKnotTable['m3_1'].augmentation_category(); A
           Aug(m3_1, GF(2))
           sage: A.cardinality()
           5
           sage: A = LegendrianKnotTable['9_13'].augmentation_category(rho=2); A
           Aug^2(9_13, GF(2))
           sage: A.poincare_polynomial(0,0) # long time (2 seconds)
           1 + 4*t
        """
        try:
            return self.__augcats__[(base_ring, rho)]
        except KeyError:
            A = AugmentationCategory_class(self, base_ring, rho)
            self.__augcats__[(base_ring,rho)] = A
            return A
    
    def augmentations(self, base_ring = GF(2), rho = 0, verbose = False):
        """
        Return a list of all augmentations of the DGA, i.e. all rho-graded
        homomorphisms self.ring() --> base_ring whose kernel contains the
        image of the differential.  Each augmentation is specified by
        a tuple of values (e(t), e(t^{-1}), e(a_0), e(a_1), ...).

        INPUT:

        - ``base_ring`` -- (default: GF(2)) A finite commutative ring
        - ``rho`` -- (default: 0) A nonnegative integer
        - ``verbose`` -- (default: False) boolean

        OUTPUT:

        A tuple of augmentations, each given as the tuple of values
        (\epsilon(x)) as x ranges over the generators of self.ring()
        starting with t and u=t^{-1}.

        EXAMPLES::

           sage: plat([2,2,2]).augmentations()
           ((1, 1, 0, 0, 1, 0, 0),
            (1, 1, 0, 1, 1, 0, 0),
            (1, 1, 1, 0, 0, 0, 0),
            (1, 1, 1, 1, 0, 0, 0),
            (1, 1, 1, 1, 1, 0, 0))
           sage: plat([2,1,4,5,1,3,5,1,2,5,4]).augmentations()
           ()
        """
        try:
            return self.__augmentations__[(base_ring, rho)]
        except KeyError:
            if self.rot() != 0 and mod(rho,2) == 0:
                self.__augmentations__[(base_ring, rho)] = ()
            else:
                R = self.__comm_ring__()
                ncab = self.__ncable__
                # Get rid of variables which must be zero for grading reasons
                zvars = list(R.gens()[:2*ncab]) # start with t_i,u_i
                vlookup = list(range(2*ncab))
                zerodict = {}
                grading_ring = IntegerModRing(rho)

                for i in range(R.ngens()-2*ncab):
                    if grading_ring(self.grading(i)) != 0:
                        zerodict[R.gen(i+2*ncab)] = 0
                    else:
                        zvars.append(R.gen(i+2*ncab))
                        vlookup.append(i+2*ncab)
                Rnew = PolynomialRing(base_ring,zvars) # only use vars with |x|=0
                tlist,ulist = Rnew.gens()[:ncab], Rnew.gens()[ncab:2*ncab]                    

                # Use Leverson's result: for Z/2-graded augs over a field, e(t)=-1
                if ncab==1 and base_ring.is_field() and rho%2==0:
                    zerodict[R.gen(0)], zerodict[R.gen(1)] = -1, -1
                    polys = [tlist[0]+1]
                polys += [tlist[i]*ulist[i]-1 for i in range(ncab)]
                
                polys += [Rnew(f.subs(zerodict)) for f in self.__abelian_differentials__()]
                nonzeropolys = list(filter(lambda p: p != 0, polys))
                
                # Find points in Spec(A_0) rather than Spec(A) -- it's much faster
                GB = Rnew.ideal(nonzeropolys).groebner_basis()
                alist = self.__augsearch__(list(GB), Rnew, base_ring)
                # Now translate points of Spec(A_0) back to points of Spec(A)
                auglist = []
                for a in alist:
                    abig = [base_ring(0)] * R.ngens()
                    for i,val in a:
                        abig[vlookup[i]] = val
                    auglist.append(tuple(abig))
                self.__augmentations__[(base_ring,rho)] = tuple(auglist)

            return self.__augmentations__[(base_ring,rho)]

    def __augsearch__(self, polys, R, coeffs, verbose=False):
        """
        Search through all augmentations of the DGA.

        INPUT:

        - ``polys`` -- a list of polynomial equations we want to solve
        - ``R`` -- the commutative ring in which these polynomials live
        - ``coeffs`` -- the underlying coefficient ring
        - ``verbose`` -- (default: False)

        OUTPUT:

        A list of augmentations.
        """
        n = R.ngens()
        genlist = list(range(n))
        
        # Sort variables by frequency so we can get rid of the most
        # common ones first
        freqlist = prod([prod(f.monomials()) for f in polys]).degrees()
        genlist.sort(key=(lambda i: -freqlist[i]))

        if verbose: import sys
        results = []

        def augsearch(glist, dlist, augmented):
            """
            Depth-first search.
            """
            if dlist == []:
                # no more equations, so remaining vars can take any values at all
                cset = FiniteEnumeratedSet(coeffs)
                if len(glist) == 0:
                    results.append(augmented)
                else:
                    for vals in cartesian_product([cset]*len(glist)):
                        results.append(augmented + list(zip(glist,vals)))
                        if verbose:
                            print(augmented)
                            sys.stdout.flush()
                return
            
            # If there are any polynomials of the form x-c,
            # make those trivial substitutions to save time
            forced, gl = [], copy(glist)
            dl = list(R.ideal(dlist).groebner_basis()) #copy(dlist)

            while True:
                subsdict = {}
                for g in copy(gl):
                    x = R.gen(g)
                    for c in coeffs:
                        if x-c in dl:
                            if x in subsdict.keys():
                                return # since some other value was already forced
                            subsdict[x] = c
                            gl.remove(g)
                            forced.append((g,c))
                newdl = [f.subs(subsdict) for f in dl]
                dl = list(filter(lambda i: i != 0, newdl))
                if len(subsdict.keys()) == 0:
                    break
                if 1 in dl:
                    return

            auglist = augmented + forced

            # Sort variables by frequency so we can get rid of the most
            # common ones first
            freqlist = R(prod([prod(f.monomials()) for f in dl])).degrees()
            gl.sort(key=(lambda i: -freqlist[i]))

            if len(gl) == 0:
                augsearch([], dl, auglist)
            else:
                x = R.gen(gl[0])
                for c in coeffs:
                    newd = [f.substitute({x:c}) for f in dl]
                    if 1 not in newd:
                        augsearch(gl[1:], list(filter(lambda d: d != 0, newd)),
                                auglist + [(gl[0],c)])

        augsearch(genlist, polys, [])
        return results

        
    def crossings(self):
        """
        Return the crossings of K.

        OUTPUT:

        A tuple of positive integers.

        EXAMPLES::

           sage: plat([2,2,2]).crossings()
           (2, 2, 2)
        """
        return self.__crossinglist__

    def differential(self, i):
        """
        Return the differential of the ith crossing of this plat.
        If i >= ncrossings, return the differential of the
        (i-ncrossings)th right cusp.

        INPUT:

        - ``i`` -- a positive integer between 1 and ncrossings+ncusps

        OUTPUT:

        An element of a free noncommutative polynomial ring.

        EXAMPLES::

           sage: K = plat([2,2,2])
           sage: K.differential(2)
           0
           sage: K.differential(3)
           t + x0 + x2 + x0*x1*x2
        """
        return self.differentials()[i]

    def differentials(self):
        """
        Return the list of differentials of the crossings and right cusps
        of this plat.

        OUTPUT:
        
        A tuple of elements of a free noncommutative polynomial ring over
        ZZ.

        EXAMPLES::

           sage: plat([2,2,2]).differentials()
           (0, 0, 0, t + x0 + x2 + x0*x1*x2, 1 - x0 - x2 - x2*x1*x0)
        """
        try:
            return self.__differential_list__
        except AttributeError:
            nc = self.ncusps()
            dlist, wlist = [],[[0]*(2*nc+1) for i in range(0,2*nc+1)]
            for i in range(nc):
                wlist[2*i+1][2*(i+1)] = 1
            ncab = self.__ncable__
            for a in range(self.ncrossings()):
                c, x = self.crossings()[a], self.ring().gen(a+2*ncab)
                gr = self.grading(a)
                dlist.append(wlist[c][c+1])
                ### Update wlist
                wlist[c][c+1] = 0
                oldwc = [wlist[i][c] for i in range(c)]
                oldwcp1 = [wlist[i][c+1] for i in range(c)]
                for i in range(1,c):
                    # w'(i,c) = w(i+1,c)+w(i,c)x, w'(i+1,c) = w(i,c)
                    wlist[i][c] = oldwcp1[i] + oldwc[i]*x
                    wlist[i][c+1] = oldwc[i]
                oldwc, oldwcp1 = copy(wlist[c]), copy(wlist[c+1])
                for j in range(c+2,2*nc+1):
                    # w'(c,j) = w(c+1,j), w'(c+1,j) = w(c,j) + x*w(c+1,j)
                    wlist[c][j] = oldwcp1[j]
                    wlist[c+1][j] = oldwc[j] - ZZ((-1)**gr)*x*oldwcp1[j]
            ### Add differentials of the cusps
            tlist = self.ring().gens()[:ncab]
            ulist = self.ring().gens()[ncab:2*ncab]
            for i in range(nc):
                if i >= ncab:
                    invisibledisk = 1
                elif self.__top_oriented_left_to_right__:
                    invisibledisk = tlist[i]
                else:
                    invisibledisk = ulist[i]
                dlist.append(invisibledisk + wlist[2*i+1][2*(i+1)])
            self.__differential_list__ = tuple(dlist)
        return self.__differential_list__

    def double_augmentation(self, e, check_augmentation = True):
        """
        Given an augmentation e of K over GF(2), return the corresponding
        augmentation of W_0(K) which augments the north and south crossings
        of each augmented crossing of K as well as the rightmost crossing
        of the clasp.

        INPUT:

        - ``e`` -- A list of nonnegative integers representing the vertices
          of K for which e(v) = 1.
        - ``check_augmentation`` -- (default: True) a boolean flag
          indicating whether we should raise a ValueError if e is not a
          valid augmentation.

        OUTPUT:

        A tuple of integers representing an augmentation of
        K.whitehead_double().

        EXAMPLES::

           sage: K = plat([2,2,2])
           sage: e = K.augmentations()[0]; K.homology(e)
           1
           sage: e2 = K.double_augmentation(e)
           (7, 8, 11, 12, 15)
           sage: K.whitehead_double().homology(e2)
           2*t + 3
        """
        R = e[0].base_ring()
        nx, nc = self.ncrossings(), self.ncusps()
        aug = [R(0)] * nc # crossings next to left cusps
        for c in e[2:nx+2]:   # skip e(t), e(u) -- they're both 1 anyway
            aug += [ R(0), c, c, R(0) ]
        return tuple([R(1),R(1)] + aug + [R(0),R(1)] + [R(0)]*(3*nc-1))

    def grading(self, n):
        """
        Return the grading of the nth crossing; if n > ncrossings, return
        the grading of the (n-ncrossings)th right cusp, counted from top
        to bottom.  This is given as an element of ZZ/(2r(K)ZZ).

        INPUT:

        - ``i`` -- an integer between 0 and ncrossings+ncusps-1

        OUTPUT:

        An element of ZZ/(2r(K)ZZ)

        EXAMPLES::

           sage: K = plat([2,2,2])
           sage: K.grading(1)
           0
           sage: K.grading(4)
           1
           sage: K.grading(4).parent()
           Integer Ring
        """
        return self.gradings()[n]

    def gradings(self):
        """
        Return a list of the gradings of this plat as elements of
        ZZ/(2r(K)ZZ).

        OUTPUT:

        A tuple of elements of ZZ/(2r(K)ZZ).

        EXAMPLES::

           sage: K = plat([2,2,2]); K.rot()
           0
           sage: K.gradings()
           (0, 0, 0, 1, 1)
           sage: K.gradings()[3].parent()
           Integer Ring
           
           sage: M = plat([2,1,4,5,1,3,5,1,2,5,4])
           sage: M.gradings()
           (0, 1, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 1, 1)
           sage: M.rot()
           -1
           sage: M.gradings()[1].parent()
           Ring of integers modulo 2
        """
        try:
            return self.__gradinglist__
        except AttributeError:
            # Strand i on the left is strand rcpos[i] on the right
            positions = list(range(2*self.ncusps()+1))
            for c in self.crossings():
                positions[c], positions[c+1] = positions[c+1], positions[c]
            rcpos = [0] + [positions.index(i)
                           for i in range(1,2*self.ncusps()+1)]
            
            # Strand i on the left has Maslov potential mu[i]
            mu = [None] * (2*self.ncusps()+1)
            #mu, cpos = [None] * (2*self.ncusps()+1), 1
            #mu[cpos] = 2*self.ncusps()
            
            # Arbitrary convention: even Maslov potiential <--> oriented from
            # left to right, odd Maslov potiential <--> right to left
            if self.__top_oriented_left_to_right__:
                cpos = rcpos.index(1)
            else:
                cpos = rcpos.index(2)
            mu[cpos] = 2*self.ncusps()
 
            while True:
                p, oldmu = rcpos[cpos], mu[cpos]
                oldpos = cpos
                if p % 2 == 0: # strand ends on the bottom of a right cusp
                    cpos = rcpos.index(p-1)
                    newmu = oldmu+1
                else:          # strand ends on the top of a right cusp
                    cpos = rcpos.index(p+1)
                    newmu = oldmu-1
                mu[cpos] = newmu

                # Move to the other strand of the new left cusp
                if cpos % 2 == 0:
                    cpos -= 1
                    newmu += 1
                else:
                    cpos += 1
                    newmu -= 1

                if mu[cpos] is not None:
                    self.__rotation_number__ = (mu[cpos]-newmu) // 2
                    break
                mu[cpos] = newmu

            if None in mu[1:]:
                raise ValueError("Plat has more than one component")

            # Use the Maslov potentials to construct the grading on crossings
            gradinglist, gring = [], Zmod(2*self.__rotation_number__)
            for c in self.crossings():
                gradinglist.append(gring(mu[c] - mu[c+1]))
                mu[c], mu[c+1] = mu[c+1], mu[c]

            # Right cusps all have grading 1
            gradinglist.extend([1] * self.ncusps())
            self.__gradinglist__ = tuple(gradinglist)
            
        return self.__gradinglist__

    def half_twist(self, ntwists = 1):
        """
        Return the 2-stranded cable of this plat constructed by adding
        n twists (n odd) to the 2-copy of this plat.  This is the satellite
        S(K,\Delta_2^n) described by Ng-Rutherford, arXiv:1206.2259.

        OUTPUT:

        A plat representing the satellite S(K,\Delta_2) of K.

        EXAMPLES::

           sage: W = plat([2,2,2]).half_twist(); W
           plat([2, 6, 4, 3, 5, 4, 4, 3, 5, 4, 4, 3, 5, 4, 1, 2, 6])
           sage: W.tb(), W.rot()
           (5, 0)
           sage: W.all_homology()
           [3]
        """
        clist = list(range(2,self.ncusps()*4, 4))
        for c in self.crossings():
            clist.extend([2*c, 2*c-1, 2*c+1, 2*c])
        clist.extend([1] * ntwists)
        clist.extend( range(2, self.ncusps()*4, 4) )

        newname = None
        if self.name() is not None:
          twiststr = '' if ntwists == 1 else '^%d'%ntwists
          newname = 'S('+self.name()+', \Delta_2'+twiststr+')'

        return plat(clist, newname)

    def homology(self, augmentation, check_augmentation = True):
        """
        Compute the linearized homology of this front with respect to the
        given graded augmentation.  The output is the *reduced* Chekanov
        polynomial p(t), i.e. the unique polynomial such that
        t + p(t) + p(t^(-1)) is the Poincare polynomial of the linearized
        homology.

        INPUT:

        - ``augmentation`` -- a list of indices of vertices satisfying
          e(v)=1
        - ``check_augmentation`` -- (default: True) a boolean flag
          determining whether or not we should raise a ValueError if the
          augmentation is not valid.

        OUTPUT:

        A polynomial in ZZ[t] representing the reduced linearized homology.

        EXAMPLES::

           sage: K = plat([2,2,2])
           sage: K.homology([1,2])
           1
           sage: K.homology([1,2]).parent()
           Univariate Polynomial Ring in t over Integer Ring
        """
        if check_augmentation and not self.is_augmentation(augmentation,0):
            raise ValueError("input is not a valid graded augmentation")
        coeffs = augmentation[0].parent() # coefficient ring, aka "base_ring"
        Rlaurent = self.__comm_ring__() # including Laurent polynomial generators
        ncab, tu_dict = self.__ncable__, {}
        for i in range(2*ncab):
            tu_dict[Rlaurent.gen(i)] = augmentation[i]
        dga_gens = Rlaurent.gens()[2*ncab:]
        n = len(dga_gens)
        R = PolynomialRing(coeffs, dga_gens)
        
        abdiffs = [R(f.subs(tu_dict)) for f in
                          self.__abelian_differentials__()]
        
        # Build lists of vertices
        gmin, gmax = min(self.gradings()), max(self.gradings())
        try:
            vdict = self.__homology_vdict__
        except AttributeError:
            vdict = {}
            vdict[gmin-1] = []
            for g in range(gmin, gmax+1):
                vdict[g] = list(filter(lambda i: self.grading(i) == g, range(n)))
            self.__homology_vdict__ = vdict

        # Twist the differential: if phi(a) = a+e(a), then d^e = phi.d.phi^{-1}
        # satisfies d^e(a) = phi(d(a-e(a))) = phi(d(a)).
        twd = []
        Phi = R.hom([R.gen(i)+augmentation[i+2*ncab] for i in range(n)])
        for d in map(Phi,abdiffs):
            # Keep only the linear terms of the twisted differential
            twd.append( R(sum([c*m if m.total_degree()==1 else 0 for c,m in d])) )

        # Build the complex from the coefficients of the linearized differentials
        ddict = {}
        for g in range(gmin,gmax+1):
            ddict[g] = matrix(coeffs, len(vdict[g-1]), len(vdict[g]),
                              lambda j,i: twd[vdict[g][i]].coefficient(R.gen(vdict[g-1][j])))

        lch_complex = ChainComplex(ddict, base_ring = coeffs, degree_of_differential=-1)
        betti = lch_complex.betti()

        S.<t> = PolynomialRing(ZZ,'t')
        positive_poly = sum([betti[d] * t**d if d>0 else 0 for d in betti.keys()])
        b_0 = self.tb() - 1 - 2*positive_poly(-1)
        reduced_poly = positive_poly - t + (b_0 // 2)
        return reduced_poly
    
    def is_augmentation(self, variable_values, rho=0):
        """
        Determine if the given list of values specify a valid
        rho-graded augmentation.

        INPUT:

        - ``variable_values`` -- the values of the augmentation, as a tuple
        - ``rho`` -- (default: 0) a nonnegative integer

        OUTPUT:

        True if the augmentation is valid and False otherwise.

        EXAMPLES::

           sage: K = plat([2,2,2])
           sage: K.is_augmentation(K.augmentations()[0])
           True
           sage: K.is_augmentation([0])
           False
        """
        if len(variable_values) != self.__comm_ring__().ngens():
            return False
        ncab = self.__ncable__
        for i in range(ncab):
            if variable_values[i]*variable_values[i+ncab] != 1: # t_i*u_i==1
                return False
        ZQ = IntegerModRing(rho)
        for i in range(self.__comm_ring__().ngens()-2*ncab):
            if variable_values[i+2*ncab] != 0 and ZQ(self.grading(i)) != 0:
                return False
        eofd = [ f(*variable_values) for f in self.__abelian_differentials__() ]
        return len(list(filter(lambda x: x != 0, eofd))) == 0

    def legendrian_mirror(self, name = None):
        """
        Return the Legendrian mirror of this knot, i.e. its image under the involution
        (x,y,z) -> (x,-y,-z), which is smoothly but not necessarily Legendrian isotopic
        to this one.

        INPUT:

        - ``name`` - (default: None) A string naming the Legendrian mirror.  If None,
        then the mirror is named 'Mirror(K)' where K is self.name().
        
        OUTPUT:

        A plat.

        EXAMPLES::

           sage: K = plat([2,4,2,3,4,4,3,2], 'Figure eight').legendrian_mirror(); K
           Mirror(Figure eight)
           sage: K.crossings()
           (4, 2, 4, 3, 2, 2, 3, 4)
        """
        Kname = self.name() if self.name() is not None else 'K'
        if name is None:
            name = 'Mirror('+Kname+')'
        mirrorc = [2*self.ncusps() - c for c in self.crossings()]
        return plat(mirrorc, name)
        
    
    def name(self):
        """
        Return the name assigned to this knot.

        OUTPUT:

        A string.

        EXAMPLES::

           sage: K = plat([2,2,2], 'trefoil')
           sage: K.name()
           'trefoil'
        """
        return self.__name__

    def ncable(self, n=2, potentials=None):
        """
        Return an n-stranded cable S(K,tw_n) of K.

        INPUT:

        - ``n`` - (default: 2) a positive integer.

        - ``potentials`` -- a list of offsets for the Maslov potentials of each strand
                            from top to bottom.  If None then the offsets are all zero.

        OUTPUT:

        A plat for the cable S(K,tw_n).

        EXAMPLES::

           sage: C = LegendrianKnotTable['m3_1'].ncable(2); C
           S(m3_1,tw_2)
           sage: C.ruling_polynomial()
           z^7 + 7*z^5 + 15*z^3 + 12*z + 4*z^-1
        """
        return _plat_cable(self, n, potentials)

    def __ncopy_differentials__(self, n):
        """
        Compute the differentials over ZZ[t^{\pm 1}] of the a^{1,n} generators
        of the n-copy using the formula of [NRSSZ, Proposition 4.14]; we store
        these with the plat in case we want to compute augmentation categories
        over lots of fields.

        INPUT:
        
        - ``n`` - a positive integer

        OUTPUT:
        
        - a tuple ((ad, xd, yd),R), where ad is a list of differentials of the a_k
          Reeb chords and xd, yd are differentials of the extra x and y
          generators corresponding to critical points on S^1.  Each differential
          is an element of the polynomial ring R.
        """
        try:
            return self.__ncopy_differential_list__[n]
        except KeyError:
            Rlaurent = self.ring() # including Laurent polynomial generators
            t,u = Rlaurent.gens()[:2]
            dga_gens = Rlaurent.gens()[2:] # Reeb chord generators
            ngen = len(dga_gens)

            # Produce a ring generated by the augmented variables
            evars = []
            for i in range(n):
                evars += ['t'+str(i), 'u'+str(i)]
                #evars += [str(dga_gens[k])+'e'+str(i) for k in range(ngen)]
                evars += ['e'+str(i)+str(dga_gens[k]) for k in range(ngen)]
            Re = PolynomialRing(ZZ, evars)

            # And a ring generated by algebra variables over the previous ring
            avars, xvars, yvars = [], [], []
            for k in range(ngen):
                akstr = str(dga_gens[k])
                avars += [akstr+'r'+str(i) for i in range(n-1)]
            for i in range(0,n-1):
                xvars.append('xr'+str(i))
                yvars.append('yr'+str(i))
            R = PolynomialRing(Re, avars + xvars + yvars)

            # Package the generators of the n-copy algebra into matrices
            # Note: anything not of the form a_k^{i,i} or a_k^{i,i+1}
            # will eventually be set to zero, so we'll save some time and
            # do it now.
            #    Doing this also lets us abelianize from the start,
            # since we'll never get products of generators of the form
            # a_{k_1}^{i_1j_1}...a_{k_m}^{i_mj_m} unless for all l, we have
            # j_l = i_{l+1} and either j_l=i_l or j_l=i_l + 1 
            MR = MatrixSpace(R,n,n)
            tm,um = MR(0), MR(0)
            for i in range(n):
                tm[i,i] = Re.gen( i*(ngen+2) )
                um[i,i] = Re.gen( i*(ngen+2) + 1)
            am = []
            for k in range(ngen):
                amk = MR(0)
                for i in range(n):
                    amk[i,i] = Re.gen( (ngen+2)*i + k+2 )
                    if i < n-1:
                        amk[i,i+1] = R.gen((n-1)*k + i)
                am.append(amk)

            xm, ym = MR(1), MR(0)
            for i in range(0,n-1):
                xm[i,i+1] = R('xr'+str(i))
                ym[i,i+1] = R('yr'+str(i))
                    
            # need to invert xm, so we'll fill in the inverse one
            # diagonal at a time
            xminv = MR(matrix(n,n,1))
            for d in range(1,n):
                for i in range(n-d): # from (0,d) to (n-d-1, n-1)
                    for j in range(i+1,i+d+1):
                        xminv[i,i+d] -= xm[i,j]*xminv[j,i+d]
            
            # The following seems to be broken for noncommutative rings,
            # so we won't use it
            #Phi = Hom(Rlaurent, MR)([tm*xm,xminv*um] + am)

            # Phi(a_k) = A_k, Phi(t) = \Delta X
            PhiDict = {t:tm*xm, u:xminv*um}
            for i in range(ngen):
                PhiDict[dga_gens[i]] = am[i]
            phidak = []

            for k in range(ngen):
                akval = MR(0)
                for m,c in Rlaurent(self.differential(k)):
                    # if the monomial is 1 there are weird substitution errors
                    # (it seems to put m.substitute(PhiDict) in a monoid rather
                    # than an algebra), so this is a hack to work around it
                    if m.is_one():
                        akval += MR(c)
                    else:
                        akval += c * Rlaurent(m).substitute(PhiDict)
                phidak.append(akval)

            dAm = []
            for k in range(ngen):
                dak = phidak[k] + ym*am[k] - ZZ((-1)**self.grading(k))*am[k]*ym
                dAm.append(dak)
            dXm = um*ym*tm*xm - xm*ym
            dYm = ym*ym

            # the only part we need is the differentials of the various c^{0,n-1}
            ndiff = tuple([m[0,n-1] for m in dAm] + [dXm[0,n-1], dYm[0,n-1]])

            # Now build a dictionary describing the linearized differentials
            # in terms of the various augmentations
            ddict, gdict = {}, {}
            Rhom = FreeAlgebra(Re, list(self.ring().gens()[2:])+['x','y'])
            for i in range(ngen):
                gdict[Rhom('a'+str(i))] = self.grading(i)+1
            gdict[Rhom('x')], gdict[Rhom('y')] = 1, 0
                
            for i in range(len(ndiff)):
                agen = Rhom.gen(i)
                for c,m in ndiff[i]:
                    # Convert a product c*a_k1^{1,2}...a_k(n-1)^{n-1,n}
                    # to c*a_k(n-1)...a_k1, or to 0 if it doesn't have this form
                    mvars = [0]*(n-1)
                    for v in m.variables():
                        a,r = str(v).split('r')
                        mvars[n-2-ZZ(r)] = Rhom(a)
                    monom = prod(mvars)

                    coeff = c
                    # compute the correct sign
                    signs = [gdict[a] for a in mvars]
                    sexp = ((n-1)*(n-2))/2 ### it's for (n-1)-fold composition
                    for p in range(n-1):
                        for q in range(p+1,n-1):
                            sexp += signs[p]*signs[q]
                        if p%2 == 1:
                            sexp += signs[p]
                    coeff *= ZZ((-1)**(sexp%2))

                    if monom in ddict.keys():
                        ddict[monom] += coeff*agen
                    else:
                        ddict[monom] = coeff*agen

            self.__ncopy_differential_list__[n] = ddict,Rhom
            return self.__ncopy_differential_list__[n]
    
    def ncrossings(self):
        """
        Return the number of crossings of this plat.

        OUTPUT:

        A positive integer.

        EXAMPLES::

           sage: plat([2,2,2]).ncrossings()
           3
        """
        return len(self.crossings())

    def ncusps(self):
        """
        Return the number of left cusps of this plat.

        OUTPUT:

        A positive integer.

        EXAMPLES::

           sage: plat([2,2,2]).ncusps()
           2
        """
        return self.__ncusps__

    def reverse(self):
        """
        Return a plat which is identical to self, but with the opposite
        orientation.

        OUTPUT:

        A plat.

        EXAMPLES::
           sage: LegendrianKnotTable['3_1'].rot()
           1
           sage: LegendrianKnotTable['3_1'].reverse().rot()
           -1
        """
        return plat(self.crossings(), self.name(),
                    top_oriented_l2r = not self.__top_oriented_left_to_right__)
    
    def ring(self):
        """
        Return a free algebra in n+2 variables over GF(2).  The differentials
        are elements of this ring modulo the relation tu=1.

        OUTPUT:

        A free algebra.

        EXAMPLES::

           sage: plat([2,2,2]).ring()
           Free Algebra on 5 generators (x0, x1, x2, x3, x4) over Finite Field of size 2
        """
        return self.__ring__

    def rot(self):
        """
        Return the rotation number of this plat for some orientation.

        OUTPUT:

        An integer.

        EXAMPLES::

           sage: plat([2,2,2]).rot()
           0
           sage: plat([2,1,4,5,1,3,5,1,2,5,4]).rot()
           -1
        """
        try:
            return self.__rotation_number__
        except AttributeError:
            _ = self.gradings()
            return self.__rotation_number__

    def ruling_polynomial(self, rho = 0, verbose = False, do_split = True):
        """
        Return the rho-graded ruling polynomial of this plat.  Note that for
        rho=0 and rho=1 this counts graded and ungraded rulings, respectively.

        INPUT:

        - ``rho'' - (default: 0) a nonnegative integer; we allow switches at
        any crossing whose grading is zero in ZZ/rho.

        OUTPUT:

        A Laurent polynomial over ZZ.

        EXAMPLES:

          sage: plat([2,2,2]).ruling_polynomial()
          z^2 + 2
          sage: LegendrianKnotTable['m7_1'].ruling_polynomial()
          z^6 + 6*z^4 + 10*z^2 + 4
          sage: LegendrianKnotTable['3_1'].ruling_polynomial()
          0
          sage: LegendrianKnotTable['3_1'].ruling_polynomial(graded=False)
          z
        """
        try:
            return self.__rp__[rho]
        except KeyError:
            from sage.rings.all import ZZ
            R.<z> = PolynomialRing(ZZ,'z')
            rz = R(0)
            ZQ = IntegerModRing(rho)
            if mod(rho,2) == 0 and self.rot() != 0:
                self.__rp__[rho] = R(0)
                return self.__rp__[rho]

            from sage.combinat.permutation import Permutation

            def _add_poly(mdict, match, poly):
              if match in mdict:
                  mdict[match] += poly
              else:
                  mdict[match] = poly

            def _is_normal(match, k):
              # True unless the pairings involving strands k,k+1 are interlaced
              tp = match[k-1]       ### the number paired with k
              bp = match[k] ### the number paired with k+1
              if k+1 < tp and tp < bp: ### k < (k+1) < k.partner < (k+1).partner
                  return False
              if tp < bp and bp < k: ### k.partner < (k+1).partner < k < k+1
                  return False
              if bp < k and k+1 < tp: ### (k+1).partner < k < k+1 < k.partner
                  return False
              return True

            nc = self.ncusps()
            rul0perm = []
            for k in range(nc):
              rul0perm += [k+k+2,k+k+1]
            rul0 = Permutation(rul0perm)
            self.__rp_counter__ = 0

            def updaterulingdict(rpdict, k, gr, ind=None):
              if verbose:
                self.__rp_counter__ += 1
                print("Crossing %d of %d, %d partial rulings; memory usage %f"%(self.__rp_counter__, self.ncrossings(), len(rpdict), get_memory_usage()))
                #print "Crossing %d of %d, %d partial rulings"%(self.__rp_counter__, self.ncrossings(), len(rpdict))
              kswap = Permutation(list(range(1,k)) + [k+1,k] + list(range(k+2,nc+nc+1)))
              rpnew = {}
              for (pm, mpoly) in rpdict.items():
                if pm[k-1] != k+1: ## "-1" because indexing starts at 0
                  ## No switching: just swap k,k+1
                  _add_poly(rpnew, kswap*pm*kswap, mpoly)
                  ## Switch if it's a normal crossing
                  if gr==0 and self.__switch_allowed__(ind) and _is_normal(pm, k):
                    _add_poly(rpnew, pm, z*mpoly)
              return rpnew

            nx = self.ncrossings()
            cglist = [(self.crossings()[i], self.grading(i)) for i in range(nx)]

            if nx < 10 or (not do_split):
              rpdict = {rul0: 1}
              for k,gr in cglist:
                rpdict = updaterulingdict(rpdict,k,ZQ(gr))
              rp_unshifted = rpdict.get(rul0, R(0))
            else:
              ## Split it in half and work from both ends
              rpd1, rpd2 = {rul0: 1}, {rul0: 1}
              fp, bp = 0, len(cglist)-1
              while fp <= bp:
                ## add a crossing to whichever dictionary is shorter
                if len(rpd1) <= len(rpd2):
                  k,gr = cglist[fp]
                  rpd1 = updaterulingdict(rpd1, k, ZQ(gr), fp)
                  fp += 1
                else:
                  k,gr = cglist[bp]
                  rpd2 = updaterulingdict(rpd2, k, ZQ(gr), bp)
                  bp -= 1
              if verbose:
                print ("Combining %d+%d partial rulings"%(len(rpd1),len(rpd2)))
              rp_unshifted = R(0)
              for (pm, mpoly) in rpd1.items():
                if pm in rpd2:
                  rp_unshifted += mpoly * rpd2[pm]

            del self.__rp_counter__

            S = LaurentPolynomialRing(ZZ,'z')
            rp = S(rp_unshifted) * S(z)**(1-nc)
            self.__rp__[rho] = rp
            return rp

    def __switch_allowed__(self, ind=None):
        """
        Auxiliary function that says whether to consider rulings with a switch at position ind.
        It can be overridden to speed up ruling polynomial computations for cables.
       """
        return True

    def tb(self):
        """
        Return the Thurston-Bennequin number of this plat.

        OUTPUT:

        An integer.

        EXAMPLES::

           sage: plat([2,2,2]).tb()
           1
           sage: plat([2,4,2,3,4,4,3,2]).tb()
           -3
        """
        try:
            return self.__tb__
        except AttributeError:
            npos = len(list(filter(lambda i: i%2 == 0, self.gradings())))
            nneg = len(self.gradings()) - npos
            self.__tb__ = npos - nneg
        return self.__tb__

    def whitehead_double(self, ntwists=0):
        """
        Return the Legendrian Whitehead double of this plat.  We add
        "ntwists" half twists next to the clasp; if ntwists is even the
        resulting plat has rotation number 0.

        INPUT:

        - ``ntwists`` - (default: 0) a nonnegative integer.

        OUTPUT:

        A plat representing the n-twisted Legendrian Whitehead double of K.

        EXAMPLES::

           sage: W = plat([2,2,2]).whitehead_double(); W
           plat([2, 6, 4, 3, 5, 4, 4, 3, 5, 4, 4, 3, 5, 4, 2, 2, 6])
           sage: W.tb(), W.rot()
           (1, 0)
           sage: W.all_homology()
           [1, 2*t + 3, t + 2]
        """
        clist = list(range(2,self.ncusps()*4, 4))
        for c in self.crossings():
            clist.extend([2*c, 2*c-1, 2*c+1, 2*c])
        clist.extend([3]*ntwists + [2])
        clist.extend( range(2, self.ncusps()*4, 4) )

        newname = None
        if self.name() is not None:
          twiststr = '' if ntwists == 0 else '_%d'%ntwists
          newname = 'Wh'+twiststr+'('+self.name()+')'

        return plat(clist, newname)

    def abelian_characteristic_algebra_nonzero(self):
        """
        Determine whether the abelianized characteristic algebra of
        this plat is nonzero.
        
        OUTPUT:
        
        Boolean.
        """
        I = ideal(self.__abelian_differentials__())
        return not (1 in I)


    def plot(self):
        """
        Print an ascii representation of the plat.

        EXAMPLES::
        
        sage: plat([2,2,2]).plot()
         __________
        /          \
        \_  _  _  _/
          \/ \/ \/
         _/\_/\_/\_
        /          \
        \__________/
        """
        self._build_plat_ascii()
        for line in self.__plat_ascii__:
            print(line)

    def _build_plat_ascii(self):
        try:
            return self.__plat_ascii__
        except AttributeError:
            self._build_left_cusps()
            for i in range(self.ncrossings()):
                crossing = self.__crossinglist__[i]
                self._build_extension()
                self._build_crossing(crossing)
            self._build_extension()
            self._build_right_cusps()
        return self.__plat_ascii__

    def _build_left_cusps(self):
        n=self.__ncusps__
        self.__plat_ascii__ = [' ']
        for i in range(n):
            self.__plat_ascii__.append('/')
            self.__plat_ascii__.append('\\')
            if i<n-1:
                self.__plat_ascii__ += [' ', ' ']

    def _build_right_cusps(self):
        for i in range(len(self.__plat_ascii__)):
            strand_char = ' '
            if i%4 == 1:
                strand_char = '\\'
            elif i%4 == 2:
                strand_char = '/'
            self.__plat_ascii__[i] = self.__plat_ascii__[i] + strand_char

    def _build_extension(self):
        for i in range(len(self.__plat_ascii__)):
            strand_char = '_' if i%2==0 else ' '
            self.__plat_ascii__[i] = self.__plat_ascii__[i] + strand_char

    def _build_crossing(self, i):
        for j in range(len(self.__plat_ascii__)):
            strand_char = '__' if j%2==0 else '  '
            if j==2*i:
                strand_char = '/\\'
            elif j==2*i-1:
                strand_char = '\\/'
            elif j==2*i-2:
                strand_char = '  '
            self.__plat_ascii__[j] = self.__plat_ascii__[j] + strand_char

###
### Augmentation categories
###

class AugmentationCategory_class:
    """
    The augmentation category of a Legendrian knot (given by a plat)
    over a finite field.

    EXAMPLES::

       sage: A = LegendrianKnotTable['m3_1'].augcat(GF(2))
       sage: A.cardinality()
       5
       sage: A = LegendrianKnotTable['m6_1'].augcat(); A
       Aug(m6_1, GF(2))
       sage: A.isomorphism_representatives() # one augmentatation per isom. class
       (0,)
       sage: A.poincare_polynomial(0,0)      # H^*hom(e_0,e_0)
       2 + t^2

       sage: # For 9_13, the graded category is empty, but the 2-graded
       sage: # category has five different isomorphism classes of objects
       sage: K = LegendrianKnotTable['9_13']
       sage: K.augcat().isomorphism_representatives()
       ()
       sage: K.augcat(rho=2).isomorphism_representatives()  # long time (4 seconds)
       (0, 2, 32, 34, 40)             
    """

    def __init__(self, K, coeffs, rho):
        if not (coeffs.is_field() and coeffs.is_finite()):
            raise ValueError("The coefficients must be a finite field")
        if not (rho in ZZ and rho >= 0):
            raise ValueError("rho must be a nonnegative integer")
        if mod(rho,2) == 1 and coeffs.characteristic() != 2:
            raise ValueError("The coefficients must have characteristic 2 unless rho is even")
        if K.rot() != 0 and rho != 0 and mod(2*K.rot(), rho) != 0:
            raise ValueError("rho must divide 2*K.rot()")
        
        self.__K__ = K
        self.__coeffs__ = coeffs
        Kname = K.name() if K.name() is not None else 'K'
        gstr = '^'+str(rho) if rho != 0 else ''
        self.__name__ = 'Aug' + gstr + '(' + Kname + ', GF(' + str(coeffs.order()) + '))'
        # free algebra for the hom spaces: generated by Reeb chords, x, y
        self.__Rhom__ = FreeAlgebra(coeffs, list(K.ring().gens()[2:])+['x','y'])
        self.__hom_complexes__, self.__hrings__, = {}, {}
        self.__composition__ = {}
        self.__n_linear__ = {}
        self.__iso_graph__ = Graph()
        self.__naut__ = {}
        self.__rho__ = rho
        if rho == 0:
            self.__ggroup__, self.__ggone__ = ZZ, ZZ(1)
        else:
            self.__ggroup__ = AdditiveAbelianGroup([self.__rho__])
            self.__ggone__ = self.__ggroup__.gen(0)       

        # gradings: we define |a^\vee| = |a|+1 for each generator a
        # we'll store them both mod 2r(K) (__gknot__) and mod rho (__glist__)
        self.__gknot__ = tuple(map(lambda g:g+1, list(map(ZZ, K.gradings())) + [0,-1]))
        self.__glist__ = tuple(map(lambda g: self.__ggone__*(g), self.__gknot__))
        
        if self.__rho__ == 0:
            self.__gmin__, self.__gmax__ = min(self.__glist__), max(self.__glist__)
            self.__grange__ = range(self.__gmin__, self.__gmax__+1)
            cxgenrange = range(self.__gmin__, self.__gmax__+2)
        else:
            self.__grange__ = self.__ggroup__
            cxgenrange = self.__ggroup__
            
        self.__gdict__ = dict([(self.__Rhom__.gen(i), self.__glist__[i])
                               for i in range(self.__Rhom__.ngens())])
        self.__cx_gens__ = {}
        for g in cxgenrange:
             self.__cx_gens__[g] = tuple(filter(lambda a: self.__gdict__[a] == g,
                                               self.__Rhom__.gens()))
        
    def __bilinearization__(self):
        """
        Return the bilinearized complex in terms of two augmentations.
        """
        try:
            return self.__bilinearization_dict__
        except AttributeError:
            diffdict, Rh = self.knot().__ncopy_differentials__(2)
            # Re == ring generated over coeffs by augmentation variables
            Re = PolynomialRing(self.coefficients(), Rh.base_ring().gens())
            Rhom = self.__Rhom__
            # Remove all the augmentation variables with nonzero grading
            Revars, nknotvars = list(Re.gens()), len(self.knot().gradings())
            for i in range(nknotvars):
                if ZZ(self.knot().grading(i)) * self.__ggone__ != 0:
                    Revars[i+2] *= 0
                    Revars[i+2+(nknotvars+2)] *= 0
            fgrading = Re.hom(Revars)

            cdict = {} # dictionary of chain complex differentials as matrices
            cgens = self.__cx_gens__

            # Build matrices representing the differentials in each grading
            ggroup, ggone = self.__ggroup__, self.__ggone__
            for g in self.__grange__:
                cdict[g] = matrix(Re, len(cgens[g+ggone]), len(cgens[g]), 0)
                for j in range(len(cgens[g])):
                    agj = Rh(cgens[g][j])
                    if agj in diffdict.keys():
                        for m,c in diffdict[agj]:
                            cg = fgrading(c) # remove e(a) where |a| is nonzero
                            if cg != 0:
                                cdict[g][cgens[g+ggone].index(Rhom(m)), j] += cg

        self.__bilinearization_dict__ = (cdict,Re)
        return self.__bilinearization_dict__

                     
    def cardinality(self):
        """
        Return the homotopy cardinality of the cohomology category, defined
        as the number of isomorphism classes of objects in which each class
        is weighted by 1/|Aut(o)| * |H_{-1}(o,o)| / |H_{-2}(o,o)| * ... 
        for a representative o.  This may not make sense if rho is 
        positive, since we need to define |Aut(o)| and the H_{-k} correctly.

        OUTPUT:
        - a rational number

        EXAMPLES::

           sage: LegendrianKnotTable['m3_1'].augcat(GF(2),rho=0).cardinality()
           5
           sage: plat([2,3]).augcat().cardinality() # a stabilized unknot
           0
        """
        try:
            return self.__cardinality__
        except AttributeError:
            self.__cardinality__ = sum([QQ(1)/self.naut(e)
                                        for e in self.isomorphism_representatives()])
        return self.__cardinality__
        
    
    def coefficients(self):
        """
        Return the coefficient field.

        OUTPUT:
        - a field

        EXAMPLES::

           sage: A = LegendrianKnotTable['1_0'].augcat(GF(5))
           sage: A.coefficients()
           Finite Field of size 5
        """
        return self.__coeffs__

    def composition(self, *augs):
        """
        Given a list of augmentations e_0, e_1, e_2, ..., e_k, return a
        function which computes the k-fold composition map
        m_k: hom(e_{k-1},e_k) \otimes ... \otimes hom(e_0,e_1) -> hom(e_0,e_k).

        INPUT:

        - ''\*augs'' -- a sequence of integers or tuples. If these are
        integers, then they refer to the augmentations self[e_0],...,self[e_k].
        It is assumed that the inputs are either all integers or all tuples.

        OUTPUT:

        - a function which takes k+1 inputs, all of which are linear
        combinations of elements of self.hom_gens(), and returns another
        such combination

        EXAMPLE:

        sage: A = plat([2,2,2]).augmentation_category(GF(3))
        sage: len(A.objects())
        10
        sage: glist = A.hom_gens(); glist
        (a0, a1, a2, a3, a4, x, y)
        sage: m2 = A.composition(0,1,2) # hom(e1,e2) x hom(e0,e1) -> hom(e0,e2)
        sage: m2(glist[0],glist[2])     # = m2(a0,a2)
        a4
        sage: m2('a3', '-y')
        a3
        """
        if augs[0] in ZZ:
            try:
                return self.__composition__[augs]
            except KeyError:
                auglist = [self[e] for e in augs]
        else:
            auglist = augs

        n = len(auglist) - 1 ## n-fold composition

        # Build a dictionary of nonzero products
        cdict, esubs, Rhom = {}, [], self.__Rhom__
        for epsilon in auglist:
            esubs += epsilon
        ddict, Rh = self.knot().__ncopy_differentials__(n+1)
        flin = Rh.base_ring().hom(esubs, self.coefficients())
        
        for a in ddict.keys():
            dda = sum([flin(c)*Rhom(m) for m,c in ddict[a]])
            if dda != 0:
                cdict[Rhom(a)] = dda

        
        def compose(*inputs):
            if len(inputs) != n:
                raise TypeError("This composition requires %s inputs (%s given)"%(n,len(inputs)))
            iprod, csum = prod(map(self.__Rhom__,inputs)), 0

            for m,c in iprod:
                mr = self.__Rhom__(m)
                if mr in cdict.keys():
                    csum += c * cdict[mr]
            return csum

        if augs[0] in ZZ:
            self.__composition__[augs] = compose
        return compose

    def __getitem__(self, n):
        """
        Return the nth object of the augmentation category with respect
        to some fixed (but arbitrary) ordering.

        INPUT:
    
        - ``n`` - An integer.
    
        OUTPUT:
    
        An augmentation of self.knot() in tuple form.

        EXAMPLES:
        
           sage: A = LegendrianKnotTable['m3_1'].augcat()
           sage: A[0]
           (1, 1, 0, 0, 1, 0, 0)
           sage: A.knot().is_augmentation(A[0])
           True
        """
        return self.objects()[n]

    def grading(self):
        """
        If this is the rho-graded augmentation category, return rho.

        OUTPUT:

        - a nonnegative integer

        EXAMPLES::
        
           sage: LegendrianKnotTable['m3_1'].augcat().grading()
           0
           sage: LegendrianKnotTable['m3_1'].augcat(rho=2).grading()
           2
        """
        return self.__rho__
    
    def hom(self, e1, e2):
        """
        Return the chain complex hom(e1,e2), where e1 and e2 are augmentations.

        INPUT:

        - ``e1``, ``e2`` - either integers or tuples.  If integers, return
        hom(self[e1],self[e2]).

        OUTPUT:

        - a ChainComplex

        EXAMPLES::

           sage: A = LegendrianKnotTable['m3_1'].augcat()
           sage: cx = A.hom(0,1); cx
           Chain complex with at most 3 nonzero terms over Finite Field of size 2
           sage: cx.betti()
           {0: 0, 1: 1, 2: 0}
        """
        if e1 in ZZ and e2 in ZZ:
            try:
                return self.__hom_complexes__[(e1,e2)]
            except KeyError:
                aug1 = self[e1]
                aug2 = self[e2]
        else:
            aug1, aug2 = e1, e2

        cd, Re = self.__bilinearization__()
        evaluate_augs = Re.hom(list(aug1) + list(aug2), self.coefficients())
        cdict = {}
        for k, dk in cd.items():
            cdict[k] = dk.apply_morphism(evaluate_augs)

        hom_complex = ChainComplex(cdict, base_ring = self.coefficients(),
                                   grading_group = self.__ggroup__,
                                   degree_of_differential=self.__ggone__)

        # Save the complex
        if e1 in ZZ and e2 in ZZ:
            self.__hom_complexes__[(e1,e2)] = hom_complex

        return hom_complex


    def hom_gens(self):
        """
        Return a tuple of variables which generate the hom spaces.

        OUTPUT:
        - a tuple of generators of a FreeAlgebra

        EXAMPLES::

           sage: gens = LegendrianKnotTable['1_0'].augcat(GF(5)).hom_gens()
           sage: gens
           (a0, x, y)
           sage: gens[0].parent()
           Free Algebra on 3 generators (a0, x, y) over Finite Field of size 5
        """
        return self.__Rhom__.gens()

    def homology_ring(self, e):
        """
        Return a FiniteDimensionalAlgebra equal to the ring H^*hom(e,e).
        Its kth basis element is called 'ek_g', where g denotes its grading.

        INPUT:
        - ``e`` -- an integer or a tuple specifying an augmentation

        OUTPUT:
        - a FiniteDimensionalAlgebra with basis elements named 'ei_j', where
        i is a nonnegative integer and j denotes the grading; if j is negative
        then we use e.g. 'n3' to denote -3

        EXAMPLES::

           sage: R = LegendrianKnotTable['m3_1'].augcat().homology_ring(0); R
           Finite-dimensional algebra of degree 3 over Finite Field of size 2
           sage: R.basis()
           [e0_0, e1_1, e2_1]
           sage: R.one()
           e0_0
           sage: e0_0, e1_1, e2_1 = R.gens()
           sage: e2_1 * (e0_0 + e1_1)
           e2_1
        """
        R,_ = self.homology_ring_with_basis(e)
        return R
        
    def homology_ring_with_basis(self, e):
        """
        Return a tuple consisting of two elements: a FiniteDimensionalAlgebra R,
        equal to the ring H^*Hom(e,e) exactly as in self.homology_ring(e);
        and a tuple of linear combinations of elements of self.hom_gens(),
        whose ith element represents the ith element of R.basis().

        INPUT:
        - ``e`` -- an integer or a tuple specifying an augmentation

        OUTPUT:
        - a tuple (R, basis), where R is a FiniteDimensionalAlgebra and basis
        is a tuple of linear combinations of elements of self.hom_gens()

        EXAMPLES::
           sage: A = LegendrianKnotTable['m3_1'].augcat(); A
           Aug(m3_1, GF(2))
           sage: R, b = A.homology_ring_with_basis(0)
           sage: R.basis()
           [e0_0, e1_1, e2_1]
           sage: b                       # e0 = y, e1 = a0+a2, e2 = a1
           (y, a0 + a2, a1)
           sage: m1 = A.composition(0,0) # the differential on A.hom(0,0)
           sage: map(m1, b)              # check: b consists of cocycles
           [0, 0, 0]
        """
        if e in ZZ and e in self.__hrings__.keys():
            return self.__hrings__[e]

        hcomplex, hgens, hbasis, hboffset = self.hom(e,e), {}, {}, {}
        H, lift, pi = {}, {}, {}
        allgens = []
        R0 = self.__Rhom__(0)

        def vector_to_ring_elt(v,g):
            cg = self.__cx_gens__[g]
            return sum([cg[i]*v[i] for i in range(len(cg))])

        def ring_elt_to_vector(x,g):
            cocycles = hcomplex.differential(g).right_kernel()
            ambient = cocycles.ambient_vector_space()
            if x==0 or g not in hgens.keys():
                return ambient(0)
            coords = [self.coefficients()(0)] * ambient.dimension()
            for m,c in x:
                coords[self.__cx_gens__[g].index(m)] = c
            return ambient(coords)
        
        # Find representatives of a basis of ker(d^g)/im(d^(g-1))
        grange = self.__grange__
        gone = self.__ggone__
            
        for g in grange:
            cg, cgl = self.__cx_gens__[g], len(self.__cx_gens__[g])
            cocycles = hcomplex.differential(g).right_kernel()
            cobdries = hcomplex.differential(g-gone).column_space()
            if cocycles.dimension() > cobdries.dimension(): # H^g != 0
                H[g], pi[g], lift[g] = cocycles.quotient_abstract(cobdries)
                hbasis[g] = [lift[g](v) for v in H[g].basis()]
                hgens[g] = [vector_to_ring_elt(v,g) for v in hbasis[g]]
                hboffset[g] = len(allgens)
                allgens += [(x,g) for x in hgens[g]]

        #unitvector = pi[0](ring_elt_to_vector(self.__Rhom__('-y'), 0))
        #unit = sum([unitvector[i] * hgens[0][i] for i in range(len(hgens[0]))])
        
        ngens = len(allgens)
        mult = self.composition(e,e,e)
        multtables, names = [], []
        for k in range(ngens):
            x,xg = allgens[k]
            gval = ZZ(str(xg).strip('()'))
            gstr = ('n' if gval<0 else '') + str(abs(gval))
            names.append('e'+str(k)+'_'+gstr)
            xmat = matrix(self.coefficients(), ngens, ngens, 0)
            # Build the matrix of right multiplication by x
            # (note: this means *rows* are entries of y*x)
            for j in range(ngens):
                y,yg = allgens[j]
                if yg+xg in hgens.keys(): # otherwise [yx] is definitely 0
                    # Convert [y*x] to a vector in terms of hbasis[g]
                    v = pi[yg+xg](ring_elt_to_vector(mult(y,x), yg+xg))
                    for i in range(len(v)):
                        if v[i] != 0:
                            xmat[j, hboffset[xg+yg]+i] = v[i]
            multtables.append(xmat)
            
        hring = FiniteDimensionalAlgebra(self.coefficients(), multtables, names)
        hrnames = tuple([d[0] for d in allgens])
        if e in ZZ:
            self.__hrings__[e] = (hring, hrnames)
        return (hring, hrnames)

    def is_isomorphic(self, e1, e2):
        """
        Determine whether two augmentations are isomorphic in the
        cohomology category H^*Aug.

        INPUT:

        - ``e1``, ``e2`` -- either both integers or both tuples.

        OUTPUT:

        - True if they are isomorphic, False otherwise.

        EXAMPLES::

           sage: A = LegendrianKnotTable['m3_1'].augcat()
           sage: A.is_isomorphic(0,1)
           False
           sage: A.is_isomorphic(3,3)
           True
        """
        if e1 in ZZ and e2 in ZZ:
            G = self.__iso_graph__
            if e1 in G.vertices() and e2 in G.vertices():
                return e2 in G.connected_component_containing_vertex(e1)
                
        # Relevant fact from [NRSSZ] section 5.3.2: an element of
        # hom(e1,e2) is an isomorphism iff it's a degree 0 cocycle
        # and its y-coefficient is a unit

        # Get a list of degree-0 cocycles which span ker(d: Hom^0 -> Hom^1)
        zero = self.__ggroup__.zero()
        bvectors = self.hom(e1,e2).differential(zero).right_kernel().basis()

        # See if one of the basis elements has nonzero y-component
        yindex = self.__cx_gens__[zero].index(self.__Rhom__('y'))
        for v in bvectors:
            if v[yindex] != 0:
                return True

        return False
    
    def isomorphism_classes(self):
        """
        Return a list of all isomorphism classes of objects in Aug.

        OUTPUT:
        - a list consisting of lists of augmentations

        EXAMPLE::
           sage: LegendrianKnotTable['m3_1'].augcat().isomorphism_classes()
           [[0], [1], [2], [3], [4]]
        """
        G = self.__iso_graph__
        for i in range(len(self.objects())):
            if i not in G.vertices():
                foundiso = -1
                for c in G.connected_components():
                    if foundiso == -1 and self.is_isomorphic(i,c[0]):
                        foundiso = c[0]
                G.add_vertex(i)
                if foundiso != -1:
                    G.add_edge(i,foundiso)
        return G.connected_components()

    def isomorphism_representatives(self):
        """
        Return a tuple of integers labeling one representative of each
        isomorphism class of object in Aug.

        OUTPUT:
        - a tuple of integers

        EXAMPLES::
           sage: A = LegendrianKnotTable['m3_1'].augcat()
           sage: A.isomorphism_representatives()
           (0, 1, 2, 3, 4)
        """
        try:
            return self.__iso_reps__
        except AttributeError:
            self.__iso_reps__ = tuple([c[0] for c in self.isomorphism_classes()])
            return self.__iso_reps__
    
    def knot(self):
        """
        Return the underlying knot.

        OUTPUT:
        - a plat

        EXAMPLES::
           sage: A = LegendrianKnotTable['5_1'].augcat()
           sage: A.knot()
           5_1
        """
        return self.__K__

    def name(self):
        """
        Return a description of this category.

        OUTPUT:
        - a string

        EXAMPLES::
           sage: A = LegendrianKnotTable['m8_21'].augcat(GF(7))
           sage: A.name()
           'Aug(m8_21, GF(7))'
        """
        return self.__name__
    
    def naut(self, e):
        """
        Return the number of automorphisms of Hom(e,e) = H^*hom(e,e),
        weighted by the negative Betti numbers
        (b_{-2}b_{-4}...) / (b_{-1}b_{-3}...)
        to account for higher homotopies; this may not make sense if
        rho is positive.

        INPUT:

        - ``e`` -- an integer or a tuple specifying an augmentation

        OUTPUT:

        - The number of automorphisms

        EXAMPLES::
        
           sage: A = LegendrianKnotTable['6_1'].augcat()
           sage: A.isomorphism_representatives()
           (0,)
           sage: A.poincare_polynomial(0)
           3 + 2*t^2
           sage: A.naut(0)
           4
        """
        if e in ZZ and e in self.__naut__.keys():
            return self.__naut__[e]

        G = self.__ggroup__
        zero, one = G.zero(), self.__ggone__
        ngen = len(self.__cx_gens__[zero])
        
        # Count the elements of ker(d_0: hom^0 --> hom^1)
        cocycles = self.hom(e,e).differential(zero).right_kernel()
        ncocycle = cocycles.cardinality()

        # Count the degree-0 cocycles with y-coefficient zero
        yindex = self.__cx_gens__[zero].index(self.__Rhom__('y'))
        proj = matrix(self.coefficients(), ngen, ngen, 0)
        proj[yindex,yindex] = 1
        nyzero = cocycles.intersection(proj.right_kernel()).cardinality()

        # Everything in the image of d_{-1}: hom^{-1} --> hom^{0}
        # has y-coefficient zero (see [NRSSZ], Lemma 5.14), so we
        # can just count all of it
        coboundaries = self.hom(e,e).differential(-1*one).column_space()
        ncobdry = coboundaries.cardinality()

        naute = (ncocycle - nyzero) // ncobdry

        q = self.coefficients().cardinality()
        for d in self.hom(e,e).betti().keys():
            dd = d if G == ZZ else d[0]
            if dd < 0:
                naute *= QQ(q) ** ((-1)**dd * self.hom(e,e).betti(d))
            
        if e in ZZ:
            self.__naut__[e] = naute
        return naute

    
    def objects(self):
        """
        Return a tuple consisting of the objects of this category, namely
        the augmentations of self.knot() over self.coefficients().

        OUTPUT:

        A tuple.

        EXAMPLES::
        
           sage: LegendrianKnotTable['1_0'].augcat().objects()
           ((1, 1, 0),)
        """
        return self.knot().augmentations(self.coefficients(),
                                         self.grading())

    def poincare_polynomial(self, e1, e2 = None):
        """
        Return the Poincare polynomial of the cohomology H^*Hom(e1,e2).

        INPUT:

        - ``e1`` -- an integer or tuple specifying an augmentation.
        - ``e2`` -- (default: None) same, but if None then we set e2=e1.

        OUTPUT:

        - A Laurent polynomial in t.

        EXAMPLES::

           sage: A = LegendrianKnotTable['m3_1'].augcat()
           sage: A.poincare_polynomial(0,0)
           1 + 2*t
           sage: A.poincare_polynomial(0,1)
           t

           sage: A = LegendrianKnotTable['7_3'].augcat()
           sage: A.isomorphism_representatives()
           (0,)
           sage: A.poincare_polynomial(0)
           2*t^-1 + 1 + 2*t^3
        """
        if e2 is None:
            e2 = e1
        bettis = self.hom(e1,e2).betti()
        S.<t> = LaurentPolynomialRing(ZZ, 't')
        poly = S(0)
        for d in bettis.keys():
            if self.__ggroup__ == ZZ:
                dd = d
            else:
                dd = 0 if self.grading() == 1 else d[0]
            poly += bettis[d] * t**dd
        return poly
        
    def __repr__(self):
        """
        Return a string describing this category.

        EXAMPLES::

           sage: LegendrianKnotTable['m8_21'].augcat(GF(7))
           Aug(m8_21, GF(7))
        """
        return self.name()
        
###
### Cables
###

class _plat_cable(plat):
    """
    Cables of plats.

    EXAMPLES::

    """
    def __init__(self, K, n, potentials = None):
        """
        Create the satellite S(K,tw_n).

        INPUT:

        - ``K`` -- a plat.

        - ``n`` -- the number of strands in the cable.

        - ``potentials`` -- a list of offsets for the Maslov potentials of each strand
                            from top to bottom.  If None then the offsets are all zero.

        EXAMPLES::

           sage: K = plat([2,2,2]) # a Legendrian trefoil (the mirror of 3_1)
           sage: C = plat_cable(K, 3)
        """
        newplat = []
        nlcusp, ncrossing, nrcusp = [], [], []
        glcusp, gcrossing, grcusp = [], [], []
        glist = []
        if potentials is None:
          potentials = [0]*n

        # build patterns for n-copies of cusps and crossings
        for i in range(n):
          ## top strand: i-1,i-2,...,0.  bottom strand: i.
          nlcusp += range(2*i, i, -1)
          glcusp += [(potentials[j]-1)-potentials[i] for j in range(i-1,-1,-1)]
          ## top strand: n-1,...,i+1. bottom strand: i.
          nrcusp += range(n+i, i+i+1, -1)
          grcusp += [potentials[j]-(potentials[i]-1) for j in range(n-1,i,-1)]
          ## top strand: n-1,...,0. bottom strand: i.
          ncrossing += range(n+i, i, -1)
          gcrossing += [potentials[j]-potentials[i] for j in range(n-1,-1,-1)]

        # add n-copies of the left cusps and crossings
        for i in range(K.ncusps()):
          newplat += [2*n*i+j for j in nlcusp]
          #glist += [-1 for _ in nlcusp]
          glist += glcusp
        # count the number of crossings involving all the n-copies of left cusps
        self.__lc_crossing_len__ = len(newplat)

        for i in range(K.ncrossings()):
          c, g = K.crossings()[i], ZZ(K.grading(i))
          newplat += [n*(c-1)+j for j in ncrossing]
          #glist += [g for _ in ncrossing]
          glist += [g+j for j in gcrossing]

        # add a full twist at the top
        for i in range(n-1):
          ## top crossing: n-i-2,...,0. bottom crossing: n-1-i.
          newplat += range(n-1, i, -1)
          #glist += [0 for _ in range(n-1, i, -1)]
          glist += [potentials[j]-potentials[n-1-i] for j in range(n-i-2,-1,-1)]
        for i in range(n-1):
          ## top crossing: i+1,...,n-1. bottom crossing: i.
          newplat += range(n-1, i, -1)
          #glist += [0 for _ in range(n-1, i, -1)]
          glist += [potentials[j]-potentials[i] for j in range(i+1,n)]

        # add the n-copies of the right cusps
        self.__rc_crossing_len__ = len(newplat)
        for i in range(K.ncusps()):
          newplat += [2*n*i+j for j in nrcusp]
          #glist += [1 for _ in nrcusp]
          glist += grcusp

        ## add gradings of right cusps
        glist += [1] * (K.ncusps() * n)

        self.__ncable__ = n
        self.__ncusps__ = K.ncusps() * n
        self.__crossinglist__ = tuple(newplat)
        ngen = len(self.__crossinglist__) + self.__ncusps__
        genstr = ','.join(['t'+str(i) for i in range(n)] +
                          ['u'+str(i) for i in range(n)] +
                          ['x'+str(j) for j in range(ngen)])
        self.__ring__ = FreeAlgebra(ZZ, n+n+ngen, genstr)
        self.__ab_ring__ = PolynomialRing(ZZ, n+n+ngen, genstr)

        self.__name__ = None
        if K.name() is not None:
          self.__name__ = 'S('+K.name()+',tw_%d)'%n

        self.__gradinglist__ = tuple(glist)
        self.__tb__ = (K.tb()+1) * n*n - n
        self.__rotation_number__ = K.rot() * n
        self.__parent_knot__ = K

        self.__augmentations__ = {}
        self.__all_homology__ = {}
        self.__ncopy_differential_list__ = {}
        self.__rp__ = {}
        self.__augcats__ = {}
        self.__top_oriented_left_to_right__ = K.__top_oriented_left_to_right__

    def augmentation_category(self, base_ring = GF(2), rho = 0):
        raise NotImplementedError("Can't do augmentation categories of cables")

    def __ncopy_differentials__(self, n):
        raise NotImplementedError("Can't compute DGA for n-copies of cables")
    
    def parent(self):
        return self.__parent_knot__

    def __switch_allowed__(self, ind=None):
        """
        Ng-Rutherford Lemma 4.6: every normal ruling of the n-cable is reduced, so there
        are no switches at crossings corresponding to the left cusps of K.
        In fact, pairs of thin strands only have switches at cusps by Lemma 3.3, so
        there are no thin pairs at all, hence no switches at n-copies of right cusps.
        """
        #return (ind is None) or ind >= self.__lc_crossing_len__
        return (ind is None) or (self.__lc_crossing_len__ <= ind and ind < self.__rc_crossing_len__)


###
###
### Some Legendrian slice knots
###
###

K9_46 = plat([2,4,6,3,5,3,4,5,5,6,2,3,4,1,2], 'm(9_46)');
K10_140_1 = plat([4,5,5,3,4,6,7,2,1,5,6,3,2,2,4,6,1,3,5,7,1,2,3,7,6,3,4], 'm(10_140)-1');
K10_140_2 = plat([2,3,1,3,4,5,2,3,4,2,1,4,3,5,1,4,2,4,3,5,4,2], 'm(10_140)-2');
K11n139 = plat([2,4,6,3,5,4,3,6,4,3,5,2,4,3,4,6,3,2,4], '11n139');
K12n582 = plat([2,1,1,4,3,2,2,5,5,4,6,3,3,2,8,7,7,6,5,4,9,8,7,6], '12n582');
K12n768 = plat([2,4,3,2,1,5,4,4,3,2,2,1,6,5,4,3,2,4,4,3,6,5,7,7,6,5,4], 'm(12n768)');
K12n838 = plat([2,4,3,3,2,1,6,5,4,3,2,7,7,6,5,6,5,4,3,2,6,5,7,3,4,7,6], '12n838');

K13n579  = plat([2,4,6,5,4,3,2,1,7,6,6,5,5,4,4,3,2,7,6,5,4,4,3,2,6,5,4], 'm(13n579)');
K13n3158 = plat([2,1,3,4,3,5,6,5,4,3,7,6,5,5,4,4,3,3,2,2,1,4,3,2,7,6,6,5,4,7,6], 'm(13n3158)');
K13n3523 = plat([2,1,3,3,2,4,6,5,2,1,5,4,3,2,8,7,6,5,4,9,8,7,6,8,7,6,5,9,8,8,7,9,9,8,7,6], '13n3523');
#K13n3523 = plat([2,1,4,6,5,5,4,3,3,2,8,7,6,5,7,6,5,4,3,8,7,6,5,5,4,7,9,8], '13n3523');
K13n4514 = plat([2,1,4,3,6,5,4,3,2,7,6,5,5,4,6,5,5,4,3,6,6,5,7,7,6,5,4], '13n4514');
K13n4659 = plat([6,5,4,3,2,1,1,3,5,4,4,3,2,2,1,5,4,3,2,7,6,6,5,4,3,2,7,6,5,5,4], 'm(13n4659)');

K14n2459 = plat([4,3,2,1,1,3,2,6,5,4,3,3,2,6,5,4,4,7,6,5,5,4,8,7,6,9,8,8], 'm(14n2459)');
K14n2601 = plat([2,4,3,6,5,4,8,7,6,5,4,3,2,2,1,9,8,8,7,6,5,4,3,2,6,5,5,4,6,9,9,8], '14n2601');
K14n8579 = plat([2,4,6,5,5,4,3,3,2,6,5,4,3,3,2,5,4,7,6], 'm(14n8579)');
K14n12406= plat([2,1,4,3,2,2,6,5,4,3,5,4,3,2,6,5,5,4,3,3,2,4,6], 'm(14n12406)');
K14n14251= plat([2,4,3,2,5,6,5,4,3,7,3,2,2,1,7,6,8,8,7,6,5,4,3,9,8,7,6,5,5,4,3,2], 'm(14n14251)');
K14n15581= plat([2,4,4,3,2,6,5,4,3,2,7,6,5,5,4,3,6,8,7,6,5,7,6,5,4,3,2,8,7,6,5,4,9,8,7,6], 'm(14n15581)');
K14n18212= plat([2,4,3,6,5,8,7,6,5,9,8,7,9,8,7,6,6,5,4,3,2,1,1,7,6,5,4,3,2,9,8,7,6,9,8,8,7,6,9,8], '14n18212');
K14n21563= plat([4,3,2,1,5,1,3,2,6,5,4,3,3,2,7,5,4,3,7,6,3,2,6,5,4,7,6], '14n21563');
K14n22150= plat([2,4,3,3,2,1,6,5,4,7,4,3,2,1,1,8,7,6,5,4,3,2,9,6,5,4,9,8,7,6,6,5,4,8,7,6], 'm(14n22150)');
K14n22789= plat([2,1,4,6,5,4,7,6,1,6,5,5,4,3,6,5,4,4,3,8,7,6,5,9,8,8,7,9,9,8,7,6,5,4,3,2], 'm(14n22789)');
K14n24246= plat([2,1,4,3,2,2,1,4,3,5,4,4,1,6,5,4,3,2,3,7,6,6,5,4,7,6,8,7,6,9,8,8], 'm(14n24246)');
K14n25967= plat([4,3,2,5,5,4,3,3,6,5,4,3,2,7,6,5,5,4,8,7,6,5,5,4,6,9,8,8], '14n25967');

# Old plats with more cusps than necessary, taken from arXiv:1107.6028
# Note -- the graded ruling polynomials of the n-cables (n \leq 3)
#         of both this 11n139 and the one above agree, likewise m(12n838)
#K11n139 = plat([4,3,3,2,1,6,5,4,3,2,2,1,1,8,7,6,5,4,3,2,5,4,4,3,9,8,8,7,6,5,4,9,9,8,7,6], '11n139');
#Km12n838 = plat([2,1,3,4,3,2,5,4,6,5,4,3,3,2,7,8,7,6,5,4,9,8,7,7,6,9,8,8,7,6,5,4,3,2,9,8,7,6,5,4], 'm(12n838)');
#K13n3523 = plat([2,4,3,6,7,8,9,10,11,12,13,14,15,16,17,5,4,3,2,4,3,7,6,5,4,6,5,9,8,7,6,11,13,15,17,8,17,16,15,14,13,12,11,10])

slice_knots = [K9_46, K10_140_1, K10_140_2, K11n139,
   K12n582,   K12n768,  K12n838,
   K13n579,   K13n3158,  K13n3523,  K13n4514,  K13n4659,
   K14n2459,  K14n2601,  K14n8579,  K14n12406, K14n14251, K14n15581,
   K14n18212, K14n21563, K14n22150, K14n22789, K14n24246, K14n25967
];

###
###
### Knot table code begins here
###
###

class __KnotTableClass__:
  def __init__(self):
    # Load the unknot
    self.__kdict__ = {'1_0': plat([], '1_0')}
    
  def __getitem__(self, val):
    """
    Returns a maximal tb representative of the knot specified by the
    string 'name'.  Acceptable names are of the form '8_19' or
    'm3_1', and all knots and their mirrors up to 9 crossings are
    included.

    INPUT:

    - ``val'' - A string of the form '8_19' or 'm3_1', or a range of such
        strings.

    OUTPUT:

    A plat of the given knot type.

    EXAMPLES::

       sage: LegendrianKnotTable['m3_1']
       m3_1
       sage: LegendrianKnotTable['7_1':'m7_1']
       [7_1, 7_2, 7_3, 7_4, 7_5, 7_6, 7_7]
       sage: LegendrianKnotTable['m3_1'].crossings()
       (2, 2, 2)
       sage: LegendrianKnotTable['m8_21'].all_homology()
       [1, t + 2]
       sage: LegendrianKnotTable['9_46'].all_homology()
       [0]
    """
    if isinstance(val,slice):
        start = (val.start == None)
        output = []
        for K in LegendrianKnotTable:
            if K.name() == val.start:
                start = True
            if K.name() == val.stop:
                return output
            if start:
                output.append(K)
        return output

    try:
        return self.__kdict__[val]
    except KeyError:
        pass

    try:    
        if val[0] == 'm':
            klist, num = __mirrors__, val[1:]
        else:
            klist, num = __knots__, val
        cr, ind = map(int, num.split('_'))
    except ValueError:
        raise ValueError("Names must be of the form '8_19' or 'm3_1'")

    if cr < 1 or cr > len(klist):
        raise ValueError("Crossing number must be between 1 and %s"%len(klist))
    if ind < 0 or ind > len(klist[cr-1]):
        raise ValueError("Index does not correspond to a %s-crossing knot"%cr)
    self.__kdict__[val] = plat(klist[cr-1][ind-1], val)
    return self.__kdict__[val]

  def __iter__(self):
    """
    Returns an iterator through the knot table.  Elements are sorted first
    by crossing number, then knots by index, then mirrors by index.
    """
    yield self['1_0']
    for i in range(2, len(__knots__)):
      for j in range(len(__knots__[i])):
        yield self['%d_%d'%(i+1,j+1)]
      for j in range(len(__mirrors__[i])):
        yield self['m%d_%d'%(i+1,j+1)]

LegendrianKnotTable = __KnotTableClass__()

__knots__ = [[], [], [[2, 1, 4, 5, 1, 3, 5, 1, 2, 5, 4]],
   [[2, 1, 4, 5, 1, 1, 3, 5, 1, 2, 5, 4]], [[4, 3, 2, 1, 6, 7, 8, 9, 
     1, 3, 5, 7, 9, 1, 2, 3, 4, 9, 8, 7, 6], [2, 3, 3, 4, 5, 1, 3, 5, 
     1, 5, 4, 3, 2]], [[2, 1, 4, 5, 3, 4, 5, 1, 4, 3, 3, 5, 1, 5, 4, 
     3, 3, 2], [2, 1, 4, 8, 7, 3, 2, 6, 7, 8, 6, 6, 5, 8, 2, 4, 7, 
     6], [2, 1, 4, 5, 6, 7, 1, 1, 3, 5, 7, 7, 6, 5, 1, 2, 5, 4]], [[6,
      5, 4, 3, 2, 1, 8, 9, 10, 11, 12, 13, 1, 3, 5, 7, 9, 11, 13, 1, 
     2, 3, 4, 5, 6, 13, 12, 11, 10, 9, 8],
    [2, 3, 3, 4, 5, 1, 3, 5, 1, 5, 4, 3, 1, 3, 2], [2, 3, 3, 4, 5, 6, 
     7, 1, 3, 5, 7, 1, 7, 6, 5, 4, 3, 3, 2], [2, 4, 6, 7, 8, 10, 11, 
     7, 6, 4, 3, 5, 4, 6, 2, 3, 9, 10, 8, 9, 10, 7, 8, 4, 5, 1, 2, 3, 
     4, 6, 8], [2, 1, 4, 6, 8, 9, 4, 5, 6, 3, 2, 4, 3, 1, 2, 6, 5, 4, 
     7, 8, 4, 6, 8],
    [6, 7, 8, 9, 4, 4, 5, 6, 7, 8, 8, 6, 4, 4, 5, 2, 3, 4, 2], [2, 1, 
     6, 8, 9, 4, 3, 2, 7, 8, 2, 4, 6, 8, 4, 5, 6, 5, 4]], [[2, 3, 4, 
     5, 1, 1, 1, 3, 5, 1, 1, 1, 5, 4, 3, 2], [4, 3, 2, 1, 6, 8, 10, 
     11, 12, 13, 6, 7, 8, 9, 10, 11, 12, 7, 6, 5, 4, 3, 2, 6, 7, 8, 7,
      6, 4, 2, 12, 10],
    [2, 1, 4, 5, 3, 2, 5, 6, 7, 8, 9, 4, 3, 1, 3, 5, 7, 9, 1, 9, 8, 7,
      6, 3, 4, 2, 3, 5, 1, 2, 5, 4], [8, 7, 6, 5, 4, 3, 10, 11, 12, 
     13, 3, 2, 1, 5, 7, 8, 9, 11, 13, 8, 7, 3, 2, 4, 3, 5, 6, 7, 1, 2,
      3, 6, 5, 3, 4, 5, 6, 11, 10, 9, 12, 11, 13, 12, 11, 10, 9, 
     8], [2, 1, 4, 5, 1, 1, 1, 3, 3, 3, 5, 5, 5, 4, 1, 2], [2, 1, 4, 
     6, 8, 4, 5, 8, 7, 6, 6, 3, 4, 2, 3, 6, 5, 4, 1, 2, 4, 6, 7, 
     8], [4, 3, 3, 2, 1, 6, 7, 7, 7, 7, 7, 3, 5, 1, 1, 2, 3, 4, 5, 
     6], [4, 3, 2, 1, 6, 7, 8, 9, 10, 11, 3, 5, 6, 9, 11, 1, 3, 7, 6, 
     5, 3, 4, 5, 11, 10, 9, 8, 7, 2, 3, 1, 2, 3, 4, 5, 6], [2, 1, 4, 
     5, 6, 7, 8, 9, 1, 1, 1, 3, 5, 7, 9, 1, 2, 3, 9, 8, 7, 6, 5, 
     4], [2, 1, 4, 6, 4, 4, 6, 6, 5, 4, 3, 2, 7, 5, 4, 2], [4, 3, 2, 
     1, 6, 8, 5, 4, 8, 7, 6, 4, 3, 2, 6, 9, 8, 6, 4, 2], [6, 7, 5, 4, 
     3, 2, 1, 1, 1, 7, 8, 9, 3, 5, 7, 9, 1, 2, 3, 4, 9, 8, 7, 6], [2, 
     1, 4, 6, 8, 9, 3, 2, 2, 7, 8, 8, 10, 10, 9, 8, 6, 4, 5, 6, 5, 4, 
     11, 10], [4, 3, 6, 7, 5, 4, 3, 2, 1, 1, 3, 4, 5, 7, 8, 9, 3, 7, 
     9, 9, 8, 7, 6, 5, 5, 4, 1, 2], [4, 6, 8, 6, 5, 8, 7, 4, 2, 1, 8, 
     6, 7, 6, 5, 4, 3, 2, 9, 8, 7, 6, 4, 2], [6, 5, 4, 3, 2, 1, 8, 7, 
     6, 5, 4, 10, 6, 5, 10, 9, 8, 7, 6, 11, 10, 8, 6, 4, 3, 2, 4, 3, 
     2, 1, 5, 4, 3, 2, 6, 4], [2, 1, 4, 3, 2, 6, 5, 7, 6, 4, 3, 2, 2, 
     8, 6, 8, 7, 9, 8, 6, 5, 4, 2, 3, 4, 3, 2, 6], [6, 5, 4, 3, 2, 1, 
     7, 8, 9, 3, 2, 1, 5, 4, 3, 3, 3, 2, 1, 7, 9, 5, 6, 7, 6, 5, 2, 3,
      1, 2, 3, 4, 9, 8, 7, 6], [2, 1, 1, 4, 3, 2, 1, 5, 4, 5, 4, 5, 4,
      5, 4, 3, 5, 4, 3, 2], [2, 1, 4, 6, 7, 1, 3, 4, 3, 1, 1, 2, 5, 6,
      4, 6], [6, 5, 8, 9, 2, 1, 4, 3, 2, 7, 5, 4, 5, 6, 7, 9, 2, 4, 5,
      6, 3, 4, 9, 8]], [[6, 5, 4, 3, 2, 1, 8, 9, 10, 11, 12, 13, 14, 
     15, 16, 17, 1, 3, 5, 7, 9, 11, 13, 15, 17, 1, 2, 3, 4, 5, 6, 17, 
     16, 15, 14, 13, 12, 11, 10, 9, 8], [2, 1, 4, 5, 1, 1, 1, 1, 1, 1,
      1, 3, 5, 1, 2, 5, 4], [4, 3, 2, 1, 6, 7, 1, 1, 3, 5, 7, 1, 1, 1,
      1, 1, 2, 3, 4, 7, 6], [6, 5, 4, 3, 2, 1, 8, 9, 1, 1, 1, 9, 7, 5,
      3, 1, 1, 1, 2, 3, 4, 9, 8, 7, 6], [8, 7, 6, 5, 4, 3, 10, 11, 12,
      13, 14, 15, 3, 5, 7, 8, 9, 11, 13, 2, 1, 8, 7, 3, 2, 1, 4, 3, 2,
      13, 12, 11, 10, 9, 5, 6, 7, 6, 5, 3, 15, 14, 13, 12, 11, 15, 14,
      13, 12, 11, 10, 9, 8, 3, 4, 5, 6], [4, 6, 8, 10, 12, 4, 5, 2, 1,
      6, 7, 12, 11, 10, 9, 8, 8, 3, 2, 4, 3, 5, 4, 6, 5, 7, 6, 1, 2, 
     3, 4, 6, 8, 13, 12, 11, 10], [2, 4, 6, 2, 3, 8, 6, 5, 4, 4, 4, 7,
      6, 9, 8, 7, 4, 5, 6, 9, 8, 1, 2, 4, 6], [2, 4, 6, 8, 10, 9, 8, 
     7, 11, 10, 9, 8, 12, 10, 9, 8, 13, 11, 10, 9, 8, 12, 11, 10, 13, 
     12, 10, 2, 3, 4, 5, 6, 6, 1, 2, 3, 4], [2, 1, 4, 6, 8, 10, 4, 3, 
     2, 1, 5, 4, 3, 2, 6, 5, 4, 3, 7, 6, 5, 4, 10, 9, 8, 8, 7, 6, 6, 
     12, 11, 10, 9, 8, 13, 12, 11, 10, 10, 13, 12], [2, 1, 4, 6, 8, 
     10, 4, 3, 2, 1, 5, 8, 7, 6, 9, 4, 3, 6, 10, 9, 8, 6, 5, 4, 4, 2, 
     8, 7, 6, 11, 10, 9, 8], [2, 3, 3, 4, 5, 6, 7, 7, 7, 3, 5, 1, 1, 
     1, 1, 2, 3, 6, 5, 5, 4], [2, 4, 6, 8, 8, 7, 6, 2, 3, 4, 5, 6, 1, 
     2, 6, 9, 8, 6, 6, 3, 4], [10, 8, 6, 4, 3, 2, 2, 1, 10, 9, 8, 7, 
     6, 5, 4, 3, 2, 11, 10, 9, 8, 4, 6, 6, 6], [2, 4, 6, 8, 10, 11, 
     12, 13, 2, 3, 4, 5, 8, 7, 6, 9, 8, 7, 10, 9, 8, 11, 10, 9, 8, 12,
      11, 10, 9, 13, 12, 10, 9, 8, 10, 8, 1, 2, 3, 4, 6], [4, 3, 3, 2,
      1, 1, 1, 3, 2, 3, 3, 6, 7, 8, 9, 10, 11, 11, 11, 11, 5, 7, 9, 
     10, 9, 8, 7, 6, 4], [4, 3, 2, 2, 2, 6, 6, 5, 4, 4, 3, 2, 6, 6, 6,
      3, 4], [6, 5, 4, 3, 8, 9, 10, 11, 12, 13, 3, 2, 1, 7, 9, 11, 10,
      12, 11, 13, 11, 1, 3, 5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 9, 10, 13, 
     12], [4, 3, 2, 2, 2, 6, 8, 8, 7, 6, 9, 8, 3, 4, 6, 6, 5, 4, 7, 6,
      2], [4, 3, 3, 2, 1, 1, 1, 6, 7, 8, 9, 9, 9, 9, 9, 7, 6, 5, 6, 7,
      7, 5, 4, 3, 2, 4, 5, 8, 7, 7, 6, 5, 4], [8, 7, 7, 6, 5, 4, 3, 3,
      2, 1, 3, 2, 10, 11, 12, 13, 13, 11, 9, 7, 5, 4, 3, 3, 13, 12, 
     11, 10, 9, 8, 7, 6, 5, 5, 4, 3, 2], [6, 5, 4, 3, 2, 1, 1, 8, 9, 
     9, 10, 11, 11, 11, 11, 7, 6, 5, 3, 6, 7, 9, 9, 8, 7, 8, 9, 1, 2, 
     3, 4, 5, 5, 6, 10, 9, 8], [4, 3, 2, 1, 5, 6, 7, 3, 7, 1, 3, 1, 3,
      5, 7, 7, 6, 5, 4, 1, 2], [4, 3, 3, 2, 1, 6, 7, 7, 7, 8, 9, 3, 2,
      1, 5, 4, 3, 2, 7, 7, 6, 5, 3, 3, 4, 9, 9, 8, 7, 9, 8, 7, 6], [2,
      4, 6, 1, 7, 8, 4, 4, 8, 7, 3, 6, 2, 2, 6, 5, 4, 8, 7, 2, 3, 6, 
     4, 3, 2], [4, 6, 3, 8, 2, 8, 1, 7, 6, 6, 5, 4, 4, 6, 9, 3, 2, 8, 
     2, 4, 6], [4, 6, 3, 7, 2, 8, 1, 9, 3, 7, 2, 6, 1, 5, 4, 3, 6, 7, 
     9, 3, 5, 7, 1, 2, 9, 8, 7, 1, 2, 7, 3, 6, 5, 1, 2, 3, 4], [2, 1, 
     4, 6, 8, 9, 6, 5, 4, 3, 2, 7, 6, 8, 7, 4, 4, 5, 6, 2, 6, 9, 8, 4,
      6], [6, 5, 4, 3, 2, 1, 8, 8, 10, 7, 11, 6, 5, 4, 8, 9, 10, 4, 3,
      2, 6, 2, 4, 8, 10], [6, 5, 4, 3, 2, 1, 1, 5, 4, 3, 3, 2, 1, 8, 
     9, 10, 11, 12, 13, 7, 6, 5, 4, 3, 2, 11, 12, 13, 9, 10, 11, 11, 
     7, 5, 3, 4, 9, 8, 7, 10, 9, 13, 12, 11, 10, 9, 8, 7, 6], [2, 4, 
     6, 8, 1, 9, 4, 3, 2, 5, 4, 3, 6, 5, 4, 1, 2, 4, 6, 7, 8, 4, 6, 8,
      8], [2, 4, 6, 8, 8, 10, 7, 9, 11, 6, 8, 10, 7, 9, 11, 8, 2, 3, 
     4, 6, 1, 4, 8, 10, 2, 5, 3, 6, 4], [4, 3, 2, 1, 6, 7, 8, 9, 1, 3,
      5, 7, 8, 1, 2, 5, 6, 7, 3, 5, 7, 9, 9, 8, 7, 7, 6, 3, 4],
    [4, 3, 2, 6, 8, 6, 5, 4, 8, 2, 4, 6, 7, 9, 3, 4, 6, 8, 6, 5, 4, 1,
      2, 7, 6],
    [4, 3, 6, 7, 3, 7, 8, 9, 2, 1, 5, 3, 5, 7, 1, 3, 4, 7, 6, 5, 9, 8,
      7, 2, 3, 9, 8, 7, 6, 1, 2, 3, 4], [2, 1, 4, 5, 1, 3, 5, 1, 3, 5,
      1, 3, 5, 5, 4, 1, 2], [2, 3, 4, 5, 6, 7, 7, 7, 1, 1, 1, 1, 5, 3,
      3, 2, 3, 6, 5, 5, 4], [2, 4, 6, 8, 8, 7, 6, 2, 3, 4, 4, 5, 6, 6,
      9, 6, 8, 5, 4, 1, 2], [2, 6, 5, 4, 8, 9, 5, 6, 6, 7, 8, 8, 8, 2,
      3, 4, 4, 4, 1, 2, 6, 5, 6, 5, 4], [6, 5, 4, 2, 1, 8, 10, 11, 8, 
     7, 6, 3, 2, 9, 8, 7, 6, 4, 5, 6, 2, 3, 4, 10, 10, 9, 8, 8, 7, 6, 
     10, 9, 8],
    [4, 3, 2, 1, 6, 5, 4, 8, 7, 9, 8, 10, 4, 6, 8, 3, 2, 4, 2, 10, 9, 
     8, 4, 5, 11, 10, 8, 7, 6, 7, 8, 5, 4],
    [4, 3, 2, 1, 6, 7, 3, 3, 1, 1, 5, 4, 3, 6, 5, 5, 7, 7, 7, 6, 1, 2,
      3, 3, 4],
    [2, 1, 6, 5, 4, 3, 2, 8, 7, 5, 4, 2, 9, 8, 6, 5, 2, 3, 9, 8, 6, 7,
      4, 3, 2, 9, 5, 6, 9, 8],
    [6, 5, 2, 1, 4, 5, 4, 3, 2, 5, 7, 4, 7, 2, 4, 7, 6], [2, 1, 4, 6, 
     8, 9, 4, 5, 6, 7, 3, 2, 1, 4, 3, 2, 6, 5, 4, 9, 7, 4, 6, 9, 
     8], [4, 3, 2, 1, 6, 7, 8, 7, 8, 6, 1, 3, 5, 6, 1, 5, 6, 5, 4, 1, 
     2], [2, 4, 6, 3, 5, 3, 4, 5, 5, 6, 2, 3, 4, 1, 2], [4, 3, 2, 1, 
     6, 7, 3, 2, 1, 5, 7, 8, 9, 2, 4, 7, 9, 1, 7, 6, 2, 3, 5, 8, 7, 9,
      8, 7, 6, 1, 2, 3, 4],
    [4, 6, 8, 9, 10, 11, 2, 1, 4, 3, 5, 7, 9, 11, 2, 4, 6, 3, 5, 7, 1,
      4, 2, 4, 6, 11, 10, 9, 8],
    [10, 9, 8, 7, 6, 2, 4, 1, 7, 8, 9, 4, 3, 5, 8, 11, 2, 4, 6, 1, 3, 
     5, 9, 2, 4, 8, 7, 11, 10, 9, 8, 4, 6]]]

__mirrors__ = [[], [], [[2, 2, 2]],
   [[2, 4, 2, 3, 4, 4, 3, 2]],
   [[2, 2, 2, 2, 2], [2, 6, 5, 4, 5, 6, 2, 3, 4, 1, 2, 4, 6]], [[4, 3,
      2, 1, 6, 5, 4, 3, 2, 8, 2, 4, 4, 8, 7, 6, 9, 8], [2, 1, 3, 4, 5,
      1, 3, 5, 1, 1, 2, 3, 5, 4], [2, 4, 6, 7, 2, 3, 4, 5, 6, 1, 2, 4,
      6, 6]], [[2, 2, 2, 2, 2, 2, 2], [2, 1, 4, 6, 8, 10, 11, 3, 2, 9,
      10, 8, 7, 6, 5, 4, 2, 10, 4, 9, 8, 7, 6], [2, 4, 6, 7, 8, 9, 5, 
     6, 7, 8, 2, 6, 1, 3, 2, 4, 6, 6, 8], [2, 1, 3, 4, 5, 1, 1, 3, 5, 
     1, 2, 3, 5, 5, 4], [2, 1, 4, 5, 6, 7, 1, 1, 1, 3, 5, 7, 7, 2, 3, 
     4, 5, 7, 6], [2, 1, 3, 4, 5, 6, 7, 7, 7, 1, 3, 5, 6, 1, 2, 3, 3, 
     3, 4], [4, 5, 6, 7, 5, 4, 3, 2, 1, 7, 6, 5, 3, 7, 5, 4, 3, 1, 1, 
     2, 3, 4, 5, 7, 6, 2, 1, 1, 2, 3, 4]], [[2, 4, 6, 8, 10, 12, 2, 3,
      4, 5, 6, 7, 12, 11, 10, 9, 8, 13, 12, 11, 10, 8, 7, 6, 1, 2, 3, 
     4],
    [2, 1, 4, 5, 5, 5, 3, 1, 1, 1, 1, 1, 5, 4, 1, 2], [2, 1, 4, 8, 7, 
     6, 7, 8, 3, 2, 4, 5, 6, 6, 6, 6, 8, 2, 3, 4], [2, 2, 2, 4, 4, 4, 
     4, 4, 3, 2, 3, 4], [4, 3, 6, 8, 10, 12, 13, 11, 10, 9, 12, 11, 
     10, 12, 11, 13, 12, 5, 4, 6, 5, 4, 3, 2, 1, 7, 6, 5, 4, 3, 2, 10,
      8, 6, 4, 3, 5, 6, 7, 4, 5, 6, 4, 5, 6, 3, 4, 10, 9, 8, 11, 10],
    [4, 3, 2, 1, 6, 7, 8, 9, 9, 1, 1, 3, 5, 7, 9, 9, 9, 8, 7, 6, 1, 2,
      3, 4], [4, 3, 2, 1, 6, 5, 4, 3, 2, 2, 4, 8, 8, 10, 10, 9, 8, 8, 
     7, 6, 11, 10, 9, 8], [2, 1, 4, 6, 7, 3, 2, 5, 6, 4, 6, 2, 6, 5, 
     4, 4, 2, 4, 7, 6], [2, 4, 6, 8, 8, 7, 6, 5, 4, 4, 3, 2, 2, 2, 9, 
     8, 7, 6, 5, 4], [4, 3, 2, 1, 6, 7, 8, 9, 10, 11, 1, 3, 5, 9, 8, 
     11, 1, 7, 8, 9, 11, 10, 7, 6, 5, 9, 8, 7, 9, 8, 7, 6, 5, 4, 1, 
     2], [4, 3, 2, 1, 6, 7, 8, 9, 1, 5, 6, 7, 8, 9, 3, 4, 5, 6, 7, 1, 
     1, 2, 3, 5, 7, 9, 3, 4, 5, 6, 9, 8],
    [2, 4, 2, 3, 6, 7, 8, 7, 6, 4, 4, 6, 5, 4, 8, 7, 6, 4, 1, 2], [2, 
     1, 4, 5, 6, 7, 3, 2, 1, 5, 6, 7, 6, 5, 5, 4, 3, 1, 7, 5, 5, 6, 7,
      6, 5, 1, 2, 3, 7, 6, 5, 4],
    [2, 4, 5, 6, 7, 8, 3, 4, 2, 2, 4, 5, 6, 8, 7, 6, 9, 8, 4, 4], [6, 
     5, 4, 3, 2, 1, 1, 1, 8, 9, 5, 4, 3, 7, 6, 5, 5, 9, 9, 9, 8, 7, 7,
      6, 1, 2, 3, 4], [2, 1, 1, 1, 4, 5, 6, 7, 3, 5, 6, 4, 5, 3, 7, 7,
      5, 1, 2, 3, 7, 6, 3, 4], [4, 3, 2, 1, 6, 7, 8, 9, 1, 3, 5, 7, 8,
      9, 6, 7, 7, 1, 3, 2, 1, 4, 3, 2, 3, 5, 4, 9, 8, 7, 7, 6],
    [4, 3, 2, 2, 1, 6, 5, 4, 4, 3, 2, 4, 8, 9, 6, 7, 8, 6, 8, 4, 5, 6,
      5, 4], [2, 4, 6, 3, 5, 7, 4, 6, 5, 3, 2, 5, 4, 6], [2, 1, 4, 3, 
     2, 2, 1, 1, 1, 6, 6, 5, 4, 3, 2, 7, 6, 5, 5, 4], [4, 3, 2, 1, 1, 
     5, 3, 2, 5, 4, 4, 4, 3, 2, 5, 4]], [[2, 2, 2, 2, 2, 2, 2, 2, 2],
    [2, 4, 6, 8, 10, 12, 14, 2, 3, 4, 5, 6, 7, 14, 13, 12, 11, 10, 9, 
     8, 8, 1, 2, 3, 4, 5, 6, 15, 14, 13, 12, 11, 10], [2, 4, 6, 8, 10,
      12, 2, 3, 4, 5, 12, 11, 10, 9, 8, 7, 6, 6, 6, 1, 2, 3, 4, 13, 
     12, 11, 10, 9, 8],
    [2, 4, 6, 8, 10, 2, 3, 4, 5, 6, 7, 10, 9, 8, 8, 8, 8, 11, 10, 1, 
     2, 3, 4, 5, 6], [2, 2, 2, 2, 4, 4, 4, 4, 4, 3, 4, 3, 2], [4, 3, 
     2, 1, 6, 7, 7, 7, 7, 7, 7, 1, 1, 5, 3, 1, 2, 3, 4, 7, 6], [6, 5, 
     4, 3, 2, 1, 8, 9, 10, 11, 1, 3, 5, 7, 9, 11, 1, 11, 11, 1, 2, 3, 
     4, 5, 6, 11, 10, 9, 8],
    [2, 1, 4, 5, 6, 7, 7, 7, 7, 7, 7, 5, 4, 3, 4, 5, 5, 6, 1, 1, 2, 3,
      3, 3, 4], [2, 1, 1, 1, 1, 1, 4, 5, 6, 7, 7, 7, 7, 7, 7, 3, 5, 6,
      2, 3, 4], [4, 3, 2, 1, 6, 7, 8, 9, 9, 1, 3, 5, 7, 9, 1, 3, 4, 9,
      9, 8, 7, 6, 5, 2, 3, 1, 2, 3, 4], [2, 1, 4, 6, 6, 8, 8, 7, 6, 5,
      4, 4, 3, 2, 9, 8, 7, 6, 4, 2], [6, 5, 4, 3, 2, 1, 8, 9, 10, 11, 
     9, 10, 11, 11, 10, 9, 1, 3, 5, 7, 9, 11, 11, 11, 10, 1, 2, 3, 4, 
     5, 6, 7, 8], [6, 5, 5, 4, 3, 2, 1, 7, 8, 9, 9, 9, 9, 9, 1, 3, 5, 
     7, 1, 2, 8, 7, 6, 5, 4], [4, 3, 6, 7, 3, 2, 1, 7, 5, 7, 1, 1, 3, 
     7, 7, 1, 2, 3, 4, 6, 7], [4, 3, 2, 1, 6, 8, 9, 6, 5, 4, 3, 2, 7, 
     8, 8, 6, 4, 3, 2, 6, 6, 5, 4, 4, 6], [8, 7, 6, 5, 7, 6, 5, 4, 3, 
     2, 1, 5, 3, 1, 1, 2, 10, 11, 12, 13, 13, 11, 9, 8, 7, 7, 13, 12, 
     11, 10, 9, 9, 8, 7, 6, 5, 4], [2, 4, 6, 4, 6, 5, 4, 4, 2, 4, 7, 
     6, 4, 3, 2, 5, 4], [6, 5, 4, 3, 2, 1, 7, 8, 9, 10, 11, 5, 3, 3, 
     4, 5, 6, 7, 1, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6, 9, 11, 11, 11, 
     10, 9, 8], [2, 1, 4, 4, 3, 2, 2, 4, 6, 6, 8, 9, 10, 11, 7, 8, 9, 
     10, 10, 6, 5, 4, 8, 7, 6], [2, 1, 4, 3, 2, 2, 2, 4, 4, 6, 6, 5, 
     4, 4, 4, 7, 6],
    [2, 1, 4, 6, 6, 8, 9, 6, 7, 8, 8, 4, 3, 2, 2, 6, 4, 5, 6, 5, 
     4], [4, 2, 2, 1, 6, 7, 8, 9, 10, 11, 12, 13, 5, 6, 7, 8, 8, 9, 
     10, 11, 12, 3, 4, 6, 4, 5, 6, 7, 8, 9, 10, 12, 9, 8, 8, 10, 2, 3,
      4, 5, 6], [2, 4, 6, 10, 1, 7, 11, 8, 4, 8, 9, 7, 6, 3, 10, 2, 6,
      4, 5, 6, 8, 10, 5, 2, 4], [6, 5, 4, 3, 2, 1, 8, 9, 10, 11, 1, 5,
      4, 3, 7, 9, 11, 3, 1, 2, 3, 3, 6, 5, 5, 4, 11, 10, 9, 8, 7, 7, 
     6], [4, 3, 2, 1, 1, 1, 1, 3, 2, 6, 5, 4, 3, 3, 7, 6, 5, 5, 4, 8, 
     7, 9, 10, 11, 11, 11, 11, 9, 10, 9, 8, 7, 6], [2, 4, 6, 7, 8, 9, 
     10, 5, 6, 2, 3, 4, 6, 7, 8, 8, 10, 1, 2, 8, 9, 10, 9, 8, 6], [6, 
     5, 4, 3, 8, 9, 10, 11, 3, 2, 1, 5, 6, 4, 5, 3, 2, 1, 7, 8, 6, 7, 
     9, 11, 2, 3, 1, 3, 5, 7, 11, 11, 10, 9, 8, 7, 6, 1, 2, 3, 4], [4,
      3, 3, 3, 3, 2, 1, 6, 7, 8, 9, 5, 4, 3, 2, 5, 4, 3, 3, 3, 9, 7, 
     5, 3, 4, 5, 6, 9, 8], [4, 3, 2, 1, 6, 5, 4, 4, 3, 2, 2, 2, 2, 6, 
     6, 6, 4, 5, 4, 3, 4, 3, 2, 7, 6], [4, 3, 3, 3, 2, 1, 6, 5, 4, 3, 
     2, 7, 8, 9, 10, 11, 7, 6, 5, 4, 3, 3, 8, 7, 7, 9, 11, 11, 11, 10,
      9, 8, 5, 3, 4, 5, 6], [4, 3, 2, 1, 6, 7, 8, 9, 3, 2, 1, 5, 7, 8,
      9, 4, 3, 9, 7, 8, 3, 5, 6, 7, 9, 6, 5, 9, 1, 2, 3, 3, 4, 5, 6, 
     9, 8], [4, 3, 2, 1, 6, 5, 4, 4, 3, 2, 2, 4, 8, 8, 10, 11, 9, 8, 
     7, 6, 10, 9, 8, 10, 9, 8, 8, 7, 6, 5, 4, 11, 10, 9, 8, 7, 6], [6,
      5, 4, 3, 2, 1, 3, 2, 1, 5, 4, 3, 3, 7, 8, 9, 10, 11, 11, 9, 7, 
     5, 6, 7, 6, 5, 1, 2, 3, 4, 11, 10, 9, 9, 8, 7, 6], [2, 1, 4, 3, 
     2, 8, 7, 6, 10, 11, 7, 8, 4, 2, 4, 5, 6, 8, 3, 9, 8, 7, 10, 9, 8,
      4, 6, 10, 2, 1, 3, 2, 4, 3, 5, 4, 6, 9, 8, 11, 10], [2, 4, 6, 8,
      9, 10, 11, 12, 13, 14, 15, 7, 8, 9, 10, 10, 11, 12, 13, 14, 14, 
     8, 12, 2, 3, 4, 5, 6, 6, 7, 8, 9, 10, 1, 2, 3, 4, 5, 6, 7, 
     8],
     ###
     ### corrected m(9_36)
     ###
     [12, 13, 6, 7, 8, 9, 10, 11, 12, 4, 5, 6, 7, 8, 2, 6, 8, 12, 2,
     1, 3, 2, 4, 6, 5, 7, 6, 8, 7, 9, 8, 10, 10],
     ###[8, 7, 6, 5, 4, 3, 2, 2, 1, 10, 9, 8, 7, 6, 5, 4, 3, 2, 4, 8,
     ### 8, 7, 6, 6, 12, 12, 11, 10, 9, 8, 13, 12], 
     [4, 3, 2, 1, 6, 7, 8,
      9, 10, 11, 9, 10, 11, 7, 8, 1, 3, 5, 9, 11, 11, 10, 7, 9, 9, 9, 
     8, 7, 6, 1, 2, 3, 4], [6, 5, 4, 3, 3, 2, 1, 5, 5, 4, 3, 2, 1, 8, 
     7, 6, 5, 4, 3, 2, 3, 9, 10, 11, 9, 10, 11, 11, 11, 8, 9, 9, 10, 
     7, 7, 8, 9, 6, 7, 8, 7, 5, 6, 5, 4], [4, 3, 2, 1, 3, 2, 1, 2, 3, 
     5, 6, 7, 8, 9, 1, 3, 5, 7, 9, 1, 3, 4, 5, 4, 3, 9, 8, 7, 6, 1, 2,
      3, 4], [4, 3, 6, 7, 3, 2, 1, 1, 3, 5, 8, 7, 9, 8, 3, 4, 5, 6, 7,
      9, 2, 3, 5, 7, 1, 2, 3, 4, 9, 8, 7, 7, 6], [2, 1, 4, 6, 8, 10, 
     11, 12, 13, 9, 10, 11, 12, 12, 8, 7, 6, 10, 9, 4, 5, 6, 8, 3, 4, 
     5, 2, 3, 4, 12, 11, 10, 9, 8, 7, 6, 1, 2, 4, 13, 12, 11, 10, 9, 
     8], [2, 1, 1, 4, 5, 3, 5, 3, 2, 4, 3, 3, 2, 4],
    [4, 3, 6, 7, 3, 2, 1, 1, 3, 5, 7, 6, 8, 4, 7, 9, 6, 8, 5, 7, 6, 8,
      7, 9, 8, 1, 2, 4, 6], [4, 3, 2, 1, 6, 7, 1, 3, 5, 7, 1, 2, 4, 7,
      4, 3, 2, 7, 6, 5, 4], [2, 4, 6, 7, 2, 3, 4, 5, 6, 4, 1, 2, 3, 6,
      5, 3, 4, 5, 7, 7, 6], [6, 5, 4, 3, 2, 1, 5, 4, 3, 2, 1, 7, 6, 7,
      2, 4, 6, 7, 1, 2, 3, 4, 5, 1, 2, 3, 4, 5, 6], [6, 7, 2, 4, 3, 3,
      4, 5, 6, 2, 4, 6, 2, 3, 4, 1, 2], [2, 1, 4, 5, 1, 3, 5, 1, 2, 5,
      4, 2, 3, 4, 1, 2, 4], [2, 1, 4, 5, 1, 3, 5, 1, 3, 5, 4, 2, 4, 3,
      2, 5, 4]]];
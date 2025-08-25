"""
Given the a Legendrian tangle in plat form together with two basepoints and a pairing of each open end, compute its commutative BLSFT algebra (see 'Maciej Wlodek. Bordered Legendrian Rational Symplectic Field Theory')
Requires lch.sage and lsft.sage, so please load those first.

AUTHORS:

- Maciej Wlodek

"""
import copy as mycopy


class BorderedPlat(BasedPlat):

    """
    A bordered plat diagram is a plat where the each of the two ends may or may not be closed up using cusps.

    EXAMPLES::

        sage: K = LegendrianKnotTable['m3_1']; B = BasedPlat(K, basepoint_bullet=(1,1), basepoint_asterisk=(2,2));
        sage: BP = BorderedPlat(B,pairing_left={1:2,3:4},pairing_right={1:2,3:4})
        sage: BP.plot()
        |____•_________|
        |              |
        |_  ___  _*_  _|
        | \/   \/   \/ |
        |_/\___/\___/\_|
        |              |
        |______________|
        |              |
        sage: BP.hamiltonian()
        t*q0*q1*q2*a1_2L*a1_2R + t*q1*q2*a1_3L*a1_2R + t*q0*q1*a1_2L*a1_3R + q0*q1*q2*a3_4L*a3_4R + t*q0*a1_2L*a1_2R + t*q2*a1_2L*a1_2R + t*q1*a1_3L*a1_3R + q0*q1*a3_4L*a2_4R + q1*q2*a2_4L*a3_4R + t*a1_3L*a1_2R + t*a1_2L*a1_3R + t*a1_4L*a1_4R + q1*a2_4L*a2_4R + q0*a3_4L*a3_4R + q2*a3_4L*a3_4R + p0*p1 + p1*p2 + p0*a2_3L + p2*a2_3R + a3_4L*a2_4R + a2_4L*a3_4R
        sage: BP.lsft_differentials()
        {t: t*b1_2L + t*b3_4L + t*b1_2R + t*b3_4R,
         u: u*b1_2L + u*b3_4L + u*b1_2R + u*b3_4R,
         q0: q0*q1*p1 + q0*q2*p2 + q0*b1_2L + q0*b1_2R + p1 + a2_3L,
         q1: q0*q1*p0 + q1*q2*p2 + q1*b1_2L + q1*b1_2R + p0 + p2,
         q2: q0*q2*p0 + q1*q2*p1 + q2*b3_4L + q2*b3_4R + p1 + a2_3R,
         p0: t*q1*q2*a1_2L*a1_2R + t*q1*a1_2L*a1_3R + q1*q2*a3_4L*a3_4R + q0*p0^2 + q1*p0*p1 + q2*p0*p2 + t*a1_2L*a1_2R + q1*a3_4L*a2_4R + p0*b1_2L + a3_4L*a3_4R + p0*b1_2R,
         p1: t*q0*q2*a1_2L*a1_2R + t*q2*a1_3L*a1_2R + t*q0*a1_2L*a1_3R + q0*q2*a3_4L*a3_4R + q0*p0*p1 + q1*p1^2 + q2*p1*p2 + t*a1_3L*a1_3R + q0*a3_4L*a2_4R + q2*a2_4L*a3_4R + p1*b1_2L + a2_4L*a2_4R + p1*b1_2R,
         p2: t*q0*q1*a1_2L*a1_2R + t*q1*a1_3L*a1_2R + q0*q1*a3_4L*a3_4R + q0*p0*p2 + q1*p1*p2 + q2*p2^2 + t*a1_2L*a1_2R + q1*a2_4L*a3_4R + p2*b3_4L + a3_4L*a3_4R + p2*b3_4R,
         a1_2L: q1*q2*a1_4L*a3_4R + q1*a1_4L*a2_4R + p0*a1_3L + a1_2L*b1_2L + a1_4L*a3_4R,
         a1_3L: q0*q1*q2*a1_4L*a3_4R + q0*q1*a1_4L*a2_4R + q0*p0*a1_3L + q1*p1*a1_3L + q2*p2*a1_3L + q0*a1_4L*a3_4R + q2*a1_4L*a3_4R + a1_2L*a2_3L + a1_4L*a2_4R + a1_3L*b1_2R,
         a1_4L: q0*p0*a1_4L + q1*p1*a1_4L + q2*p2*a1_4L + a1_2L*a2_4L + a1_3L*a3_4L + a1_4L*b3_4L + a1_4L*b1_2R,
         a2_3L: t*q0*q1*q2*a1_3L*a1_2R + t*q0*q1*a1_3L*a1_3R + q0*q1*q2*a2_4L*a3_4R + t*q0*a1_3L*a1_2R + t*q2*a1_3L*a1_2R + q0*q1*a2_4L*a2_4R + q0*p0*a2_3L + q1*p1*a2_3L + q2*p2*a2_3L + t*a1_3L*a1_3R + q0*a2_4L*a3_4R + q2*a2_4L*a3_4R + a2_3L*b1_2L + a2_4L*a2_4R + a2_3L*b1_2R,
         a2_4L: t*q0*q1*q2*a1_4L*a1_2R + t*q0*q1*a1_4L*a1_3R + t*q0*a1_4L*a1_2R + t*q2*a1_4L*a1_2R + q0*p0*a2_4L + q1*p1*a2_4L + q2*p2*a2_4L + t*a1_4L*a1_3R + a2_3L*a3_4L + a2_4L*b1_2L + a2_4L*b3_4L + a2_4L*b1_2R,
         a3_4L: t*q1*q2*a1_4L*a1_2R + t*q1*a1_4L*a1_3R + t*a1_4L*a1_2R + p0*a2_4L + a3_4L*b3_4L,
         b1_2L: t*q1*q2*a1_3L*a1_2R + t*q1*a1_3L*a1_3R + q1*q2*a2_4L*a3_4R + t*a1_3L*a1_2R + t*a1_4L*a1_4R + q1*a2_4L*a2_4R + p0*a2_3L + b1_2L^2 + a2_4L*a3_4R,
         b3_4L: t*q1*q2*a1_3L*a1_2R + t*q1*a1_3L*a1_3R + q1*q2*a2_4L*a3_4R + t*a1_3L*a1_2R + t*a1_4L*a1_4R + q1*a2_4L*a2_4R + p0*a2_3L + b3_4L^2 + a2_4L*a3_4R,
         a1_2R: q0*q1*a3_4L*a1_4R + q1*a2_4L*a1_4R + p2*a1_3R + a3_4L*a1_4R + a1_2R*b1_2R,
         a1_3R: q0*q1*q2*a3_4L*a1_4R + q1*q2*a2_4L*a1_4R + q0*p0*a1_3R + q1*p1*a1_3R + q2*p2*a1_3R + q0*a3_4L*a1_4R + q2*a3_4L*a1_4R + b3_4L*a1_3R + a2_4L*a1_4R + a1_2R*a2_3R + a1_3R*b1_2R + a1_3R*b3_4R,
         a1_4R: q0*p0*a1_4R + q1*p1*a1_4R + q2*p2*a1_4R + b3_4L*a1_4R + a1_2R*a2_4R + a1_3R*a3_4R + a1_4R*b1_2R,
         a2_3R: t*q0*q1*q2*a1_2L*a1_3R + t*q1*q2*a1_3L*a1_3R + q0*q1*q2*a3_4L*a2_4R + t*q0*a1_2L*a1_3R + t*q2*a1_2L*a1_3R + q1*q2*a2_4L*a2_4R + t*a1_3L*a1_3R + q0*p0*a2_3R + q1*p1*a2_3R + q2*p2*a2_3R + q0*a3_4L*a2_4R + q2*a3_4L*a2_4R + b3_4L*a2_3R + a2_4L*a2_4R + a2_3R*b3_4R,
         a2_4R: t*q0*q1*q2*a1_2L*a1_4R + t*q1*q2*a1_3L*a1_4R + t*q0*a1_2L*a1_4R + t*q2*a1_2L*a1_4R + t*a1_3L*a1_4R + q0*p0*a2_4R + q1*p1*a2_4R + q2*p2*a2_4R + b3_4L*a2_4R + a2_3R*a3_4R,
         a3_4R: t*q0*q1*a1_2L*a1_4R + t*q1*a1_3L*a1_4R + t*a1_2L*a1_4R + p2*a2_4R + a3_4R*b3_4R,
         b1_2R: t*q0*q1*a1_2L*a1_3R + t*q1*a1_3L*a1_3R + q0*q1*a3_4L*a2_4R + t*a1_2L*a1_3R + t*a1_4L*a1_4R + q1*a2_4L*a2_4R + p2*a2_3R + a3_4L*a2_4R + b1_2R^2,
         b3_4R: t*q0*q1*a1_2L*a1_3R + t*q1*a1_3L*a1_3R + q0*q1*a3_4L*a2_4R + t*a1_2L*a1_3R + t*a1_4L*a1_4R + q1*a2_4L*a2_4R + p2*a2_3R + a3_4L*a2_4R + b3_4R^2}
    """

    def __init__(self, based_plat=None, pairing_left=None, pairing_right=None):
        """
        Create a bordered plat from a based plat. pairing_left and pairing_right specify how to connect the left and right endpoints together to form the knot. (see Section 5.1 of 'Bordered Legendrian Rational Symplectic Field Theory')
        Warning: the program currently does not check that the pairings are valid (i.e. form a single-component knot). Leave either pairing as None to close off that end with cusps instead of a dividing line.
        

        INPUT:

        - ``based_plat`` -- an object of class BasedPlat from which to construct the BorderedPlat.

        - ``pairing_left`` -- a dictionary with integer keys and integer values, assumed to all be in the range [1,2*num_strands]. 
          This specifies how to connect the left endpoints together. Leave as None to instead close up the left endpoints with left cusps.

        - ``pairing_right`` -- a dictionary with integer keys and integer values, assumed to all be in the range [1,2*num_strands]. 
          This specifies how to connect the right endpoints together. Leave as None to instead close up the right endpoints with right cusps.


        EXAMPLES::

            sage: K = LegendrianKnotTable['m3_1']; B = BasedPlat(K, basepoint_bullet=(1,1), basepoint_asterisk=(2,2));
            sage: BorderedPlat(B,pairing_left={1:2,3:4},pairing_right={1:2,3:4}).plot()
            |____•_________|
            |              |
            |_  ___  _*_  _|
            | \/   \/   \/ |
            |_/\___/\___/\_|
            |              |
            |______________|
            |              |
            sage: BorderedPlat(B,pairing_left=None,pairing_right={1:2,3:4}).plot()
             ____•_________|
            /              |
            \_  ___  _*_  _|
              \/   \/   \/ |
             _/\___/\___/\_|
            /              |
            \______________|
                           |
            sage: BorderedPlat(B,pairing_left={1:2,3:4},pairing_right=None).plot()
            |____•_________ 
            |              \
            |_  ___  _*_  _/
            | \/   \/   \/  
            |_/\___/\___/\_ 
            |              \
            |______________/
            |       
            sage: BorderedPlat(B,pairing_left=None,pairing_right=None).plot()
             ____•_________ 
            /              \
            \_  ___  _*_  _/
              \/   \/   \/  
             _/\___/\___/\_ 
            /              \
            \______________/
        """


        super().__init__(crossinglist=based_plat.crossings(), basepoint_bullet=based_plat.basepoint_bullet, basepoint_asterisk=based_plat.basepoint_asterisk, swap_basepoints=based_plat.swap_basepoints, name=based_plat.name)
        self.pairing_left=pairing_left
        self.pairing_right=pairing_right

        c = len(self.__crossinglist__)
        nc = self.__ncusps__
        n = c + nc if pairing_right is None else c

        self.n_cr = n

        lsft_varnames = ['t','u'] 
        lsft_varnames += ['q0'] if n==1 else ['q'+str(i) for i in range(n)]
        lsft_varnames += ['p0'] if n==1 else ['p'+str(i) for i in range(n)]
        if pairing_left is not None:
            for i in range(1,2*nc+1):
                for j in range(i+1,2*nc+1):
                    lsft_varnames += ['a'+str(i)+'_'+str(j)+'L']
            for i in pairing_left:
                j = pairing_left[i]
                if j>i and j <= 2*self.__ncusps__:
                    lsft_varnames += ['b'+str(i)+'_'+str(j)+'L']
        if pairing_right is not None:
            for i in range(1,2*nc+1):
                for j in range(i+1,2*nc+1):
                    lsft_varnames += ['a'+str(i)+'_'+str(j)+'R']
            for i in pairing_right:
                j = pairing_right[i]
                if j>i and j <= 2*self.__ncusps__:
                    lsft_varnames += ['b'+str(i)+'_'+str(j)+'R']

        

        self.__lsft_ring__ = PolynomialRing(GF(2), len(lsft_varnames),lsft_varnames)


    def plot(self):
        """
        Print an ascii representation of the bordered plat.

        EXAMPLES::

            sage: K = LegendrianKnotTable['m4_1']; B = BasedPlat(K, basepoint_bullet=(3,2), basepoint_asterisk=(7,4))
            sage: BP = BorderedPlat(B, pairing_right={1:2,3:4,5:6}); BP.plot()
             _____________________________|
            /                             |
            \_  ____  _•_______________  _|
              \/    \/                 \/ |
             _/\____/\___  _______  ___/\_|
            /            \/       \/      |
            \____  ______/\_  _  _/\_*____|
                 \/         \/ \/         |
             ____/\_________/\_/\_________|
            /                             |
            \_____________________________|
                                          |
        """

        self._build_based_plat_ascii()
        for line in self.__plat_ascii__:
            print(line)


    def hamiltonian(self):
        """
        Return the hamiltonian of this bordered plat.

        OUTPUT:
        
        An element of self.lsft_ring().

        EXAMPLES::

            sage: K = LegendrianKnotTable['m4_1']; B = BasedPlat(K, basepoint_bullet=(3,2), basepoint_asterisk=(7,4))
            sage: BP = BorderedPlat(B, pairing_right={1:2,3:4,5:6}); BP.hamiltonian()
            t*q0*q1*q6*p7 + t*q4*q5*p2*p7 + t*q3*q6*p2*p7 + u*q0*q2*q7*a1_2R + q0*q1*q2*p3 + q4*q5*p3*p6 + t*q0*q5*p7 + q0*q4*q5*a1_2R + q0*q3*q6*a1_2R + u*q0*q2*a1_3R + t*q0*q1*a3_4R + t*q3*p2*a3_4R + t*q4*p2*a3_5R + q1*q4*q7*a3_6R + q1*q4*q6*a4_6R + q1*q4*q5*a5_6R + q0*q2*p4 + q3*p1*p4 + q5*p1*p6 + t*p2*p7 + u*q7*a1_2R + q0*q3*a1_4R + q0*q4*a1_5R + q1*q4*a2_6R + q7*p6*a3_4R + t*q0*a3_5R + q7*p5*a3_5R + q3*q7*a3_6R + q4*p3*a4_5R + q6*p5*a4_5R + q3*q6*a4_6R + q3*q5*a5_6R + p0*p2 + q1*p3 + p4*p5 + p3*p6 + q0*a1_2R + u*a1_3R + p7*a2_3R + p6*a2_4R + p5*a2_5R + q3*a2_6R + p1*a4_5R + q1*a5_6R + p4 + a4_6R
        """

        try:
            return self.__hamiltonian__
        except AttributeError:
            nc = self.ncusps()
            n = self.n_cr
            alpha = [[0]*(2*nc+1) for i in range(0,2*nc+1)]
            hamiltonian = 0
            if self.pairing_left is None:
                for i in range(nc):
                    alpha[2*i+1][2*(i+1)] = 1
            else:
                for i in range(1,2*nc+1):
                    for j in range(i+1,2*nc+1):
                        alpha[i][j] = self.lsft_ring()('a'+str(i)+'_'+str(j)+'L')
            bx = self.basepoint_bullet[0]
            bz = self.basepoint_bullet[1]
            t = self.lsft_ring()('t')
            u = self.lsft_ring()('u')
            # """%TODO t,u
            if bx == 0:
                if self.pairing_left is None:
                    if bz % 2 == 0:
                        alpha[bz-1][bz] *= u
                    else:
                        alpha[bz][bz+1] *= t
                else:
                    for i in range(1,bx):
                        alpha[i][bx] *= u
                    for i in range(bx+1,2*nc+1):
                        alpha[bx][i] *= t
            # """
            for a in range(self.ncrossings()):
                c, q, p = self.crossings()[a], self.lsft_ring().gen(a+2), self.lsft_ring().gen(a+2+n)
                #for row in alpha:
                #    print(row)
                #print('-------------')
                old_alpha = mycopy.deepcopy(alpha)
                for i in range(1, c):
                    alpha[i][c] = old_alpha[i][c+1] + old_alpha[i][c]*q
                    alpha[i][c+1] = old_alpha[i][c]

                hamiltonian = hamiltonian + old_alpha[c][c+1]*p
                alpha[c][c+1] = p
                for j in range(c+2,2*nc+1):
                    alpha[c][j] = old_alpha[c+1][j]
                    alpha[c+1][j] = old_alpha[c][j] + q*old_alpha[c+1][j]

                # """TODO t,u
                if bx == a+1:
                    for i in range(1,bz):
                        alpha[i][bz] = alpha[i][bz]*u
                    for i in range(bz+1, 2*nc+1):
                        alpha[bz][i] = alpha[bz][i]*t
                # """
            
            #for row in alpha:
            #    print(row)
            #print('-------------')

            if self.pairing_right is None:
                for i in range(nc):
                    p = self.lsft_ring().gen(i+2+self.ncrossings()+n)
                    hamiltonian = hamiltonian + p + p*alpha[2*i+1][2*i+2]
            else:
                for i in range(1,2*nc+1):
                    for j in range(i+1,2*nc+1):
                        aR = self._alpha(i,j,'R')
                        hamiltonian += alpha[i][j]*aR

            self.__alpha__ = alpha
            self.__hamiltonian__ = hamiltonian
            return self.__hamiltonian__
    
    def bracket(self, x, y):
        """
        Return the SFT bracket of x and y.

        INPUT:

        - ``x`` -- an element of self.lsft_ring()

        - ``y`` -- an element of self.lsft_ring()

        OUTPUT:
        
        An element of self.lsft_ring().

        EXAMPLES::

            sage: K = LegendrianKnotTable['m4_1']; B = BasedPlat(K, basepoint_bullet=(3,2), basepoint_asterisk=(7,4))
            sage: BP = BorderedPlat(B, pairing_right={1:2,3:4,5:6});
            sage: x=BP.lsft_ring()('a1_2R')
            sage: y=BP.lsft_ring()('a2_4R')
            sage: BP.bracket(x,y)
            a1_4R
        """

        n = self.n_cr
        nc = self.__ncusps__

        x_grad = x.gradient()
        # print(x_grad)
        # print(len(x_grad))
        # y_grad = self.dual(y).gradient()
        y_grad = y.gradient()
        # print(y_grad)
        # print(len(y_grad))
        bracket_xy = 0
        # for i in range(2, 2*n+2):
        #     bracket_xy += x_grad[i] * self.dual(y_grad[i])
        for i in range(2,n+2):
            bracket_xy += x_grad[i] * y_grad[i+n] + x_grad[i+n] * y_grad[i]
        k = 2*n+2

        if self.pairing_left is not None:
            (bracket_xy, k) = self._alpha_beta_bracket(x_grad,y_grad,k,nc,bracket_xy,'L')
        if self.pairing_right is not None:
            (bracket_xy, k) = self._alpha_beta_bracket(x_grad,y_grad,k,nc,bracket_xy,'R')

        return bracket_xy


    def lsft_differentials(self):
        """
        Return the lsft differentials of the generators of the lsft algebra associated to this bordered plat.

        OUTPUT:
        
        A dictionary whose keys are the generators of self.lsft_ring() and whose values are polynomials in self.lsft_ring().

        EXAMPLES::

            sage: K = LegendrianKnotTable['m4_1']; B = BasedPlat(K, basepoint_bullet=(3,2), basepoint_asterisk=(7,4))
            sage: BP = BorderedPlat(B, pairing_right={1:2,3:4,5:6});
            sage: BP.lsft_differentials()
            {t: t*b1_2R + t*b3_4R + t*b5_6R,
             u: u*b1_2R + u*b3_4R + u*b5_6R,
             q0: q0*q2*p2 + q0*q3*p3 + q0*q4*p4 + q0*q5*p5 + q0*q6*p6 + q0*q7*p7 + q0*b1_2R + p2,
             q1: q1*q4*p4 + q1*q5*p5 + q3*p4 + q5*p6 + q1*b5_6R + a4_5R,
             q2: t*q4*q5*p7 + t*q3*q6*p7 + q0*q2*p0 + q2*q3*p3 + q2*q4*p4 + q2*q5*p5 + q2*q6*p6 + q2*q7*p7 + t*q3*a3_4R + t*q4*a3_5R + t*p7 + q2*b1_2R + p0,
             q3: q0*q1*q2 + q0*q3*p0 + q2*q3*p2 + q3^2*p3 + q3*q4*p4 + q3*q5*p5 + q4*q5*p6 + q4*a4_5R + q3*b5_6R + q1 + p6,
             q4: q0*q4*p0 + q1*q4*p1 + q2*q4*p2 + q3*q4*p3 + q4^2*p4 + q0*q2 + q3*p1 + p5 + 1,
             q5: q0*q5*p0 + q1*q5*p1 + q2*q5*p2 + q3*q5*p3 + q5^2*p5 + q7*a3_5R + q6*a4_5R + p4 + a2_5R,
             q6: q0*q6*p0 + q2*q6*p2 + q4*q5*p3 + q6^2*p6 + q5*p1 + q7*a3_4R + q6*b5_6R + p3 + a2_4R,
             q7: t*q0*q1*q6 + t*q4*q5*p2 + t*q3*q6*p2 + t*q0*q5 + q0*q7*p0 + q2*q7*p2 + t*p2 + q7*b1_2R + a2_3R,
             p0: t*q1*q6*p7 + u*q2*q7*a1_2R + q0*p0^2 + q2*p0*p2 + q1*q2*p3 + q3*p0*p3 + q4*p0*p4 + q5*p0*p5 + q6*p0*p6 + t*q5*p7 + q7*p0*p7 + q4*q5*a1_2R + q3*q6*a1_2R + u*q2*a1_3R + t*q1*a3_4R + q2*p4 + q3*a1_4R + q4*a1_5R + t*a3_5R + p0*b1_2R + a1_2R,
             p1: t*q0*q6*p7 + q1*p1^2 + q0*q2*p3 + q4*p1*p4 + q5*p1*p5 + t*q0*a3_4R + q4*q7*a3_6R + q4*q6*a4_6R + q4*q5*a5_6R + q4*a2_6R + p1*b5_6R + p3 + a5_6R,
             p2: u*q0*q7*a1_2R + q0*p0*p2 + q2*p2^2 + q0*q1*p3 + q3*p2*p3 + q4*p2*p4 + q5*p2*p5 + q6*p2*p6 + q7*p2*p7 + u*q0*a1_3R + q0*p4 + p2*b1_2R,
             p3: t*q6*p2*p7 + q0*p0*p3 + q2*p2*p3 + q4*p3*p4 + q5*p3*p5 + q0*q6*a1_2R + t*p2*a3_4R + p1*p4 + q0*a1_4R + q7*a3_6R + q6*a4_6R + q5*a5_6R + p3*b5_6R + a2_6R,
             p4: t*q5*p2*p7 + q0*p0*p4 + q1*p1*p4 + q2*p2*p4 + q3*p3*p4 + q5*p3*p6 + q0*q5*a1_2R + t*p2*a3_5R + q1*q7*a3_6R + q1*q6*a4_6R + q1*q5*a5_6R + q0*a1_5R + q1*a2_6R + p3*a4_5R,
             p5: t*q4*p2*p7 + q0*p0*p5 + q1*p1*p5 + q2*p2*p5 + q3*p3*p5 + q4*p3*p6 + t*q0*p7 + q0*q4*a1_2R + q1*q4*a5_6R + p1*p6 + q3*a5_6R,
             p6: t*q0*q1*p7 + t*q3*p2*p7 + q0*p0*p6 + q2*p2*p6 + q0*q3*a1_2R + q1*q4*a4_6R + p5*a4_5R + q3*a4_6R + p6*b5_6R,
             p7: u*q0*q2*a1_2R + q0*p0*p7 + q2*p2*p7 + q7*p7^2 + q1*q4*a3_6R + u*a1_2R + p6*a3_4R + p5*a3_5R + q3*a3_6R + p7*b1_2R,
             a1_2R: q1*q4*a1_6R + p7*a1_3R + p6*a1_4R + p5*a1_5R + q3*a1_6R + a1_2R*b1_2R,
             a1_3R: t*q0*q1*a1_4R + t*q3*p2*a1_4R + t*q4*p2*a1_5R + q1*q4*q7*a1_6R + q0*p0*a1_3R + q2*p2*a1_3R + q7*p7*a1_3R + q7*p6*a1_4R + t*q0*a1_5R + q7*p5*a1_5R + q3*q7*a1_6R + a1_2R*a2_3R,
             a1_4R: q1*q4*q6*a1_6R + q0*p0*a1_4R + q2*p2*a1_4R + q7*p7*a1_4R + q4*p3*a1_5R + q6*p5*a1_5R + q3*q6*a1_6R + p1*a1_5R + a1_2R*a2_4R + a1_3R*a3_4R + a1_4R*b3_4R + a1_6R,
             a1_5R: q1*q4*q5*a1_6R + q0*p0*a1_5R + q1*p1*a1_5R + q2*p2*a1_5R + q3*p3*a1_5R + q6*p6*a1_5R + q7*p7*a1_5R + q3*q5*a1_6R + q1*a1_6R + a1_2R*a2_5R + a1_3R*a3_5R + a1_4R*a4_5R + a1_5R*b1_2R,
             a1_6R: q0*p0*a1_6R + q1*p1*a1_6R + q2*p2*a1_6R + q3*p3*a1_6R + q6*p6*a1_6R + q7*p7*a1_6R + a1_2R*a2_6R + a1_3R*a3_6R + a1_4R*a4_6R + a1_5R*a5_6R + a1_6R*b1_2R + a1_6R*b5_6R,
             a2_3R: u*q0*q2*q7*a1_3R + q0*q4*q5*a1_3R + q0*q3*q6*a1_3R + t*q0*q1*a2_4R + t*q3*p2*a2_4R + t*q4*p2*a2_5R + q1*q4*q7*a2_6R + u*q7*a1_3R + q0*p0*a2_3R + q2*p2*a2_3R + q7*p7*a2_3R + q7*p6*a2_4R + t*q0*a2_5R + q7*p5*a2_5R + q3*q7*a2_6R + q0*a1_3R + a2_3R*b1_2R,
             a2_4R: u*q0*q2*q7*a1_4R + q0*q4*q5*a1_4R + q0*q3*q6*a1_4R + q1*q4*q6*a2_6R + u*q7*a1_4R + q0*p0*a2_4R + q2*p2*a2_4R + q7*p7*a2_4R + q4*p3*a2_5R + q6*p5*a2_5R + q3*q6*a2_6R + q0*a1_4R + p1*a2_5R + a2_3R*a3_4R + a2_4R*b1_2R + a2_4R*b3_4R + a2_6R,
             a2_5R: u*q0*q2*q7*a1_5R + q0*q4*q5*a1_5R + q0*q3*q6*a1_5R + q1*q4*q5*a2_6R + u*q7*a1_5R + q0*p0*a2_5R + q1*p1*a2_5R + q2*p2*a2_5R + q3*p3*a2_5R + q6*p6*a2_5R + q7*p7*a2_5R + q3*q5*a2_6R + q0*a1_5R + q1*a2_6R + a2_3R*a3_5R + a2_4R*a4_5R,
             a2_6R: u*q0*q2*q7*a1_6R + q0*q4*q5*a1_6R + q0*q3*q6*a1_6R + u*q7*a1_6R + q0*p0*a2_6R + q1*p1*a2_6R + q2*p2*a2_6R + q3*p3*a2_6R + q6*p6*a2_6R + q7*p7*a2_6R + q0*a1_6R + a2_3R*a3_6R + a2_4R*a4_6R + a2_5R*a5_6R + a2_6R*b5_6R,
             a3_4R: u*q0*q2*a1_4R + q1*q4*q6*a3_6R + q4*p3*a3_5R + q6*p5*a3_5R + q3*q6*a3_6R + u*a1_4R + p7*a2_4R + p1*a3_5R + a3_4R*b3_4R + a3_6R,
             a3_5R: u*q0*q2*a1_5R + q1*q4*q5*a3_6R + q1*p1*a3_5R + q3*p3*a3_5R + q6*p6*a3_5R + q3*q5*a3_6R + u*a1_5R + p7*a2_5R + q1*a3_6R + a3_4R*a4_5R + a3_5R*b1_2R,
             a3_6R: u*q0*q2*a1_6R + q1*p1*a3_6R + q3*p3*a3_6R + q6*p6*a3_6R + u*a1_6R + p7*a2_6R + a3_4R*a4_6R + a3_5R*a5_6R + a3_6R*b1_2R + a3_6R*b5_6R,
             a4_5R: t*q0*q1*a3_5R + t*q3*p2*a3_5R + q1*q4*q5*a4_6R + q0*q3*a1_5R + q7*p6*a3_5R + q1*p1*a4_5R + q3*p3*a4_5R + q6*p6*a4_5R + q3*q5*a4_6R + p6*a2_5R + q1*a4_6R + a4_5R*b1_2R + a4_5R*b3_4R,
             a4_6R: t*q0*q1*a3_6R + t*q3*p2*a3_6R + q0*q3*a1_6R + q7*p6*a3_6R + q1*p1*a4_6R + q3*p3*a4_6R + q6*p6*a4_6R + p6*a2_6R + a4_5R*a5_6R + a4_6R*b1_2R + a4_6R*b3_4R + a4_6R*b5_6R,
             a5_6R: t*q4*p2*a3_6R + q0*q4*a1_6R + t*q0*a3_6R + q7*p5*a3_6R + q4*p3*a4_6R + q6*p5*a4_6R + p5*a2_6R + p1*a4_6R + a5_6R*b5_6R,
             b1_2R: u*q0*q2*a1_3R + q0*q3*a1_4R + q0*q4*a1_5R + q1*q4*a2_6R + u*a1_3R + p7*a2_3R + p6*a2_4R + p5*a2_5R + q3*a2_6R + b1_2R^2,
             b3_4R: u*q0*q2*a1_3R + t*q4*p2*a3_5R + q1*q4*q7*a3_6R + q1*q4*q6*a4_6R + q0*q3*a1_4R + t*q0*a3_5R + q7*p5*a3_5R + q3*q7*a3_6R + q4*p3*a4_5R + q6*p5*a4_5R + q3*q6*a4_6R + u*a1_3R + p7*a2_3R + p6*a2_4R + p1*a4_5R + b3_4R^2 + a4_6R,
             b5_6R: t*q4*p2*a3_5R + q1*q4*q7*a3_6R + q1*q4*q6*a4_6R + q0*q4*a1_5R + q1*q4*a2_6R + t*q0*a3_5R + q7*p5*a3_5R + q3*q7*a3_6R + q4*p3*a4_5R + q6*p5*a4_5R + q3*q6*a4_6R + p5*a2_5R + q3*a2_6R + p1*a4_5R + b5_6R^2 + a4_6R}
        """

        try:
            return self.__lsft_differentials__
        except AttributeError:
            self.__lsft_differentials__ = dict()
            hamiltonian = self.hamiltonian()
            string_differentials = self.string_differentials()
            n = self.n_cr
            nc = self.__ncusps__
            # for i in range(n):
            #     p_i = self.lsft_ring().gen(i+2+n)
            #     q_i = self.lsft_ring().gen(i+2)
            #     self.__lsft_differentials__[p_i] = self.bracket(hamiltonian, p_i) + string_differentials[p_i]
            #     self.__lsft_differentials__[q_i] = self.bracket(hamiltonian, q_i) + string_differentials[q_i]
            for gen in self.lsft_ring().gens():
                self.__lsft_differentials__[gen] = self.bracket(hamiltonian, gen) + string_differentials[gen]

            return self.__lsft_differentials__

    def string_differentials(self):
        """
        Return the string differentials of the generators of the lsft algebra associated to this bordered plat.

        OUTPUT:
        
        A dictionary whose keys are the generators of self.lsft_ring() and whose values are polynomials in self.lsft_ring().

        EXAMPLES::

            sage: K = LegendrianKnotTable['m4_1']; B = BasedPlat(K, basepoint_bullet=(3,2), basepoint_asterisk=(7,4))
            sage: BP = BorderedPlat(B, pairing_right={1:2,3:4,5:6});
            sage: BP.string_differentials()
            {q0: q0*q2*p2 + q0*q3*p3 + q0*q4*p4 + q0*q5*p5 + q0*q6*p6 + q0*q7*p7 + q0*b1_2R,
             p0: q0*p0^2 + q2*p0*p2 + q3*p0*p3 + q4*p0*p4 + q5*p0*p5 + q6*p0*p6 + q7*p0*p7 + p0*b1_2R,
             q1: q1*q4*p4 + q1*q5*p5 + q1*b5_6R,
             p1: q1*p1^2 + q4*p1*p4 + q5*p1*p5 + p1*b5_6R,
             q2: q0*q2*p0 + q2*q3*p3 + q2*q4*p4 + q2*q5*p5 + q2*q6*p6 + q2*q7*p7 + q2*b1_2R,
             p2: q0*p0*p2 + q2*p2^2 + q3*p2*p3 + q4*p2*p4 + q5*p2*p5 + q6*p2*p6 + q7*p2*p7 + p2*b1_2R,
             q3: q0*q3*p0 + q2*q3*p2 + q3^2*p3 + q3*q4*p4 + q3*q5*p5 + q3*b5_6R,
             p3: q0*p0*p3 + q2*p2*p3 + q4*p3*p4 + q5*p3*p5 + p3*b5_6R,
             q4: q0*q4*p0 + q1*q4*p1 + q2*q4*p2 + q3*q4*p3 + q4^2*p4,
             p4: q0*p0*p4 + q1*p1*p4 + q2*p2*p4 + q3*p3*p4,
             q5: q0*q5*p0 + q1*q5*p1 + q2*q5*p2 + q3*q5*p3 + q5^2*p5,
             p5: q0*p0*p5 + q1*p1*p5 + q2*p2*p5 + q3*p3*p5,
             q6: q0*q6*p0 + q2*q6*p2 + q6^2*p6 + q6*b5_6R,
             p6: q0*p0*p6 + q2*p2*p6 + p6*b5_6R,
             q7: q0*q7*p0 + q2*q7*p2 + q7*b1_2R,
             p7: q0*p0*p7 + q2*p2*p7 + q7*p7^2 + p7*b1_2R,
             t: t*b1_2R + t*b3_4R + t*b5_6R,
             u: u*b1_2R + u*b3_4R + u*b5_6R,
             a1_2R: a1_2R*b1_2R,
             a1_3R: q0*p0*a1_3R + q2*p2*a1_3R + q7*p7*a1_3R + a1_2R*a2_3R,
             a1_4R: q0*p0*a1_4R + q2*p2*a1_4R + q7*p7*a1_4R + a1_2R*a2_4R + a1_3R*a3_4R + a1_4R*b3_4R,
             a1_5R: q0*p0*a1_5R + q1*p1*a1_5R + q2*p2*a1_5R + q3*p3*a1_5R + q6*p6*a1_5R + q7*p7*a1_5R + a1_2R*a2_5R + a1_3R*a3_5R + a1_4R*a4_5R + a1_5R*b1_2R,
             a1_6R: q0*p0*a1_6R + q1*p1*a1_6R + q2*p2*a1_6R + q3*p3*a1_6R + q6*p6*a1_6R + q7*p7*a1_6R + a1_2R*a2_6R + a1_3R*a3_6R + a1_4R*a4_6R + a1_5R*a5_6R + a1_6R*b1_2R + a1_6R*b5_6R,
             a2_3R: q0*p0*a2_3R + q2*p2*a2_3R + q7*p7*a2_3R + a2_3R*b1_2R,
             a2_4R: q0*p0*a2_4R + q2*p2*a2_4R + q7*p7*a2_4R + a2_3R*a3_4R + a2_4R*b1_2R + a2_4R*b3_4R,
             a2_5R: q0*p0*a2_5R + q1*p1*a2_5R + q2*p2*a2_5R + q3*p3*a2_5R + q6*p6*a2_5R + q7*p7*a2_5R + a2_3R*a3_5R + a2_4R*a4_5R,
             a2_6R: q0*p0*a2_6R + q1*p1*a2_6R + q2*p2*a2_6R + q3*p3*a2_6R + q6*p6*a2_6R + q7*p7*a2_6R + a2_3R*a3_6R + a2_4R*a4_6R + a2_5R*a5_6R + a2_6R*b5_6R,
             a3_4R: a3_4R*b3_4R,
             a3_5R: q1*p1*a3_5R + q3*p3*a3_5R + q6*p6*a3_5R + a3_4R*a4_5R + a3_5R*b1_2R,
             a3_6R: q1*p1*a3_6R + q3*p3*a3_6R + q6*p6*a3_6R + a3_4R*a4_6R + a3_5R*a5_6R + a3_6R*b1_2R + a3_6R*b5_6R,
             a4_5R: q1*p1*a4_5R + q3*p3*a4_5R + q6*p6*a4_5R + a4_5R*b1_2R + a4_5R*b3_4R,
             a4_6R: q1*p1*a4_6R + q3*p3*a4_6R + q6*p6*a4_6R + a4_5R*a5_6R + a4_6R*b1_2R + a4_6R*b3_4R + a4_6R*b5_6R,
             a5_6R: a5_6R*b5_6R,
             b1_2R: b1_2R^2,
             b3_4R: b3_4R^2,
             b5_6R: b5_6R^2}
        """

        try:
            return self.__string_differentials__
        except AttributeError:
            n = self.n_cr
            nc = self.__ncusps__

            betas = dict()
            holomorphic_corners = dict()

            betas_L = dict()
            betas_R = dict()

            for i in range(n):
                q = self.lsft_ring().gen(2+i)
                betas[q] = 0
                holomorphic_corners[q] = False

            for i in range(1,2*self.__ncusps__+1):
                betas_L[i] = 0
                betas_R[i] = 0

            start_x = self.basepoint_bullet[0]
            start_z = self.basepoint_bullet[1]

            #first go right
            self._traverse(betas, betas_L, betas_R, holomorphic_corners, start_x, start_z, True)
    
            #then go left
            self._traverse(betas, betas_L, betas_R, holomorphic_corners, start_x, start_z, False)

            # print(betas_L)
            # print(betas_R)

            self.__string_differentials__ = dict()
            for i in range(n):
                q = self.lsft_ring().gen(2+i)
                p = self.lsft_ring().gen(2+n+i)
                self.__string_differentials__[q] = betas[q]*q if holomorphic_corners[q] else betas[q]*q+q*p*q
                self.__string_differentials__[p] = betas[q]*p+p*q*p if holomorphic_corners[q] else betas[q]*p

            
            t = self.lsft_ring()('t')
            u = self.lsft_ring()('u')
            self.__string_differentials__[t] = 0
            self.__string_differentials__[u] = 0

            if self.pairing_left is not None:
                for i in range(1,2*nc+1):
                    for j in range(i+1, 2*nc+1):
                        alpha = self._alpha(i,j,'L')
                        self.__string_differentials__[alpha] = (betas_L[i] + betas_L[j])*alpha
                        for m in range(i+1,j):
                            self.__string_differentials__[alpha] += self._alpha(i,m,'L') * self._alpha(m,j,'L')
                for i in self.pairing_left:
                    j = self.pairing_left[i]
                    if j>i and j <= 2*self.__ncusps__:
                        beta = self._beta(i,j,'L')
                        self.__string_differentials__[beta] = beta*beta
                        self.__string_differentials__[t] += beta
                        self.__string_differentials__[u] += beta
            if self.pairing_right is not None:
                for i in range(1,2*nc+1):
                    for j in range(i+1, 2*nc+1):
                        alpha = self._alpha(i,j,'R')
                        self.__string_differentials__[alpha] = (betas_R[i] + betas_R[j])*alpha
                        for m in range(i+1,j):
                            self.__string_differentials__[alpha] += self._alpha(i,m,'R') * self._alpha(m,j,'R')
                for i in self.pairing_right:
                    j = self.pairing_right[i]
                    if j>i and j <= 2*self.__ncusps__:
                        beta = self._beta(i,j,'R')
                        self.__string_differentials__[beta] = beta*beta
                        self.__string_differentials__[t] += beta
                        self.__string_differentials__[u] += beta

            self.__string_differentials__[t] *= t
            self.__string_differentials__[u] *= u

            return self.__string_differentials__


    def _traverse(self, betas, betas_L, betas_R, holomorphic_corners, x, z, right=True):
        target_x = self.basepoint_asterisk[0]
        target_z = self.basepoint_asterisk[1]
        target_reached = (x == target_x) and (z == target_z) and (self.swap_basepoints == right)
        beta = 0
        while not target_reached:
            (beta, x, z, right, target_reached) = self._step(betas, betas_L, betas_R, holomorphic_corners, beta, x, z, target_x, target_z, right)
            #print('beta = {}, x = {}, z = {}, right = {}, target_reached = {}'.format(beta, x, z, right, target_reached))

    #return tuple (new_beta, new_x, new_z, new_right, target_reached)
    def _step(self, betas, betas_L, betas_R, holomorphic_corners, beta, x, z, target_x, target_z, right=True):
        #right cusp reached - update cusp and turn back.
        if x == self.ncrossings() and right:
            if self.pairing_right is None:
                q = self.lsft_ring().gen(2+self.ncrossings()+(z-1)//2)
                p = self.lsft_ring().gen(2+self.ncrossings()+self.ncusps()+self.ncrossings()+(z-1)//2)
                betas[q] = q*p
                holomorphic_corners[q] = True
                new_z = self._new_z_cusp(z)
                return (beta, x, new_z, False, False)
            else:
                new_z = self._new_z_pairing(self.pairing_right, z)
                betas_R[z] = beta
                beta += self._beta(z, new_z, 'R')
                betas_R[new_z] = beta
                return (beta, x, new_z, False, False)
        #left cusp reached
        elif x == 0 and not right:
            if self.pairing_left is None:
                new_z = self._new_z_cusp(z)
                target_reached = (x==target_x) and ((z==target_z) or (new_z==target_z))
                return (beta, x, new_z, True, target_reached)
            else:
                new_z = self._new_z_pairing(self.pairing_left, z)
                if x == target_x and z == target_z:
                    return (beta, x, new_z, True, True)
                target_reached = (x==target_x) and ((z==target_z) or (new_z==target_z))
                betas_L[z] = beta
                beta += self._beta(z,new_z,'L')
                if target_reached:
                    return (beta, x, new_z, True, target_reached)
                betas_L[new_z] = beta
                return (beta, x, new_z, True, target_reached)

        else:
            if right:
                crossing = self.__crossinglist__[x]
                new_x = x+1
                new_z = self._new_z_crossing(z, crossing)
                target_reached = (new_x == target_x) and (new_z == target_z)
                if z!=new_z:
                    q = self.lsft_ring().gen(2+x)
                    p = self.lsft_ring().gen(2+self.n_cr+x)
                    betas[q] = betas[q] + beta
                    beta = beta + q*p
                return (beta, new_x, new_z, True, target_reached)
            else:
                target_reached = (x == target_x) and (z == target_z)
                if target_reached:
                    return (beta, x, z, False, True)
                else:
                    new_x = x-1
                    crossing = self.__crossinglist__[new_x]
                    new_z = self._new_z_crossing(z, crossing)
                    if z!=new_z:
                        q = self.lsft_ring().gen(2+new_x)
                        p = self.lsft_ring().gen(2+self.n_cr+new_x)
                        betas[q] = betas[q] + beta
                        beta = beta+q*p
                        holomorphic_corners[q] = not holomorphic_corners[q]
                    return (beta, new_x, new_z, False, False)

    def _new_z_crossing(self, z, crossing):
        return z+1 if z==crossing else (z-1 if z==crossing+1 else z)

    def _new_z_cusp(self, z):
        return (z-1) if z%2==0 else (z+1)

    def _new_z_pairing(self, pairing, z):
        for (i,j) in pairing.items():
            if z == i:
                return j
            elif z == j:
                return i
        return -1

    def _build_left_cusps(self):
        cusps = self.pairing_left is None

        n=self.__ncusps__
        self.__plat_ascii__ = [' '] if cusps else ['|']
        for i in range(n):
            self.__plat_ascii__.append('/' if cusps else '|')
            self.__plat_ascii__.append('\\' if cusps else '|')
            if i<n-1:
                self.__plat_ascii__ += [' ', ' '] if cusps else ['|', '|']
        self.__plat_ascii__.append(' ' if cusps else '|')

    def _build_right_cusps(self):
        cusps = self.pairing_right is None

        for i in range(len(self.__plat_ascii__)):
            if not cusps:
                strand_char = '|'
            else:
                strand_char = ' '
                if i%4 == 1:
                    strand_char = '\\'
                elif i%4 == 2:
                    strand_char = '/'
            self.__plat_ascii__[i] = self.__plat_ascii__[i] + strand_char

    def _alpha_beta_bracket(self, x_grad,y_grad,k,nc, bracket_xy, L='L'):
        x_alpha = [[0]*(2*nc+1) for i in range(0,2*nc+1)]
        y_alpha = [[0]*(2*nc+1) for i in range(0,2*nc+1)]
        x_alpha_sum = [0]*(2*nc+1)
        y_alpha_sum = [0]*(2*nc+1)

        x_beta = [[0]*(2*nc+1) for i in range(0,2*nc+1)]
        y_beta = [[0]*(2*nc+1) for i in range(0,2*nc+1)]

        beta_dict = self.pairing_left if L=='L' else self.pairing_right

        c = k
        for i in range(1,2*nc+1):
            for j in range(i+1,2*nc+1):
                # print('{},{},{}'.format(i,j,c))
                x_alpha[i][j] = x_grad[c]
                y_alpha[i][j] = y_grad[c]
                x_alpha_sum[j] += x_grad[c] * self._alpha(i,j,L)
                x_alpha_sum[i] += x_grad[c] * self._alpha(i,j,L)
                y_alpha_sum[j] += y_grad[c] * self._alpha(i,j,L)
                y_alpha_sum[i] += y_grad[c] * self._alpha(i,j,L)
                c += 1
        for i in beta_dict:
            j = beta_dict[i]
            if j>i and j <= 2*self.__ncusps__:
                x_beta[i][j] = x_grad[c]
                y_beta[i][j] = y_grad[c]
                c += 1


        for i in range(1,2*nc+1):
            for j in range(i+1,2*nc+1):
                for m in range(i+1,j):
                    bracket_xy += (x_alpha[i][m]*y_alpha[m][j] + x_alpha[m][j]*y_alpha[i][m]) * self._alpha(i,j,L)
        for i in beta_dict:
            j = beta_dict[i]
            if j>i and j <= 2*self.__ncusps__:
                bracket_xy += (x_alpha_sum[j] + x_alpha_sum[i])*y_beta[i][j]
                bracket_xy += (y_alpha_sum[j] + y_alpha_sum[i])*x_beta[i][j]
        return (bracket_xy, c)

    def _alpha(self, i, j, L='L'):
        if i>j:
            temp=j
            j=i
            i=temp
        return self.__lsft_ring__('a'+str(i)+'_'+str(j)+L)

    def _beta(self, i, j, L='L'):
        if i>j:
            temp=j
            j=i
            i=temp
        return self.__lsft_ring__('b'+str(i)+'_'+str(j)+L)
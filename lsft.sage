"""
Given the a Legendrian front in plat form together with two basepoints, compute its commutative LSFT algebra (see 'Lenhard Ng. Rational Symplectic Field Theory for Legendrian Knots. Invent. Math., 182(3):451–512, 2010.')
Requires lch.sage, so please load that first.

AUTHORS:

- Maciej Wlodek

"""


import copy as mycopy

class BasedPlat(plat):

    """
    A based plat diagram is a plat with two specified basepoints, represented by a bullet and an asterisk.

    EXAMPLES::

        sage: K = LegendrianKnotTable['m3_1']; B = BasedPlat(K, basepoint_bullet=(1,1), basepoint_asterisk=(2,2)); B.plot()
         ____•_________
        /              \
        \_  ___  _*_  _/
          \/   \/   \/
         _/\___/\___/\_
        /              \
        \______________/
        sage: B.hamiltonian()
        t*q0*q1*q2*p3 + q0*q1*q2*p4 + t*q0*p3 + t*q2*p3 + p0*p1 + p1*p2 + q0*p4 + q2*p4 + p3 + p4
        sage: B.lsft_differentials()
        {p0: t*q1*q2*p3 + q0*p0^2 + q1*p0*p1 + q2*p0*p2 + q1*q2*p4 + t*p3 + p4,
        q0: q0*q1*p1 + q0*q2*p2 + p1,
        p1: t*q0*q2*p3 + q0*p0*p1 + q1*p1^2 + q2*p1*p2 + q0*q2*p4,
        q1: q0*q1*p0 + q1*q2*p2 + p0 + p2,
        p2: t*q0*q1*p3 + q0*p0*p2 + q1*p1*p2 + q2*p2^2 + q0*q1*p4 + t*p3 + p4,
        q2: q0*q2*p0 + q1*q2*p1 + p1,
        p3: 0,
        q3: t*q0*q1*q2 + q3^2*p3 + t*q0 + t*q2 + 1,
        p4: 0,
        q4: q0*q1*q2 + q4^2*p4 + q0 + q2 + 1,
        t: 0,
        u: 0}
    """

    def __init__(self, plat=None, basepoint_bullet=(0,1), basepoint_asterisk=(0,1), swap_basepoints=False, crossinglist=None, name=None, **extra_args):
        """
        Create a based plat, either from a plat (default) or from a crossing list. 
        The two basepoints are each specified by a pair of integers (x,z), where x 
        is the index of the crossing after which the basepoint should be placed, while 
        z is the number of the strand on which the basepoint should be placed.
        If x is 0, the basepoint is placed before the first crossing.
        x must be in the range [0,num_crossings], while z must be in the range [1,num_strands] = [1,2*num_cusps].

        INPUT:

        - ``plat`` -- an object of class Plat from which to construct the BasedPlat.

        - ``basepoint_bullet -- a 2-tuple of integers describing the position of the 'bullet' basepoint. 
          The first coordinate specifies the position along the x-axis and must be in [0, num_crossings].
          The second coordinate specifies the position along the z-axis and must be in [1, num_strands].

        - ``basepoint_asterisk -- a 2-tuple of integers describing the position of the 'asterisk' basepoint. 
          The first coordinate specifies the position along the x-axis and must be in [0, num_crossings].
          The second coordinate specifies the position along the z-axis and must be in [1, num_strands].
          
        - ``swap_basepoints`` -- a boolean which is applicable only if basepoint_bullet = basepoint_asterisk.
          The default behavior in such a scenario is to put the 'bullet' basepoint immediately to the right of 
          the 'asterisk' basepoint. However, if swap_basepoints is set to True, the 'bullet' basepoint is placed
          immediately to the left of the 'asterisk' basepoint.

        - ``crossinglist`` -- an optional list of positive integers describing the crossings from left to right. 
          Used in lieu of plat only if plat is not specified.

        - ``name`` -- an optional string such as '3_1'.

        EXAMPLES::

            sage: K=LegendrianKnotTable['m3_1']
            sage: BasedPlat(K).plot()
             _*_•__________ 
            /              \
            \_____  _  _  _/
                  \/ \/ \/  
             _____/\_/\_/\_ 
            /              \
            \______________/
            sage: BasedPlat(K,swap_basepoints=True).plot()
             _•_*__________ 
            /              \
            \_____  _  _  _/
                  \/ \/ \/  
             _____/\_/\_/\_ 
            /              \
            \______________/
            sage: BasedPlat(crossinglist=[2,2,2], basepoint_bullet=(0,4), basepoint_asterisk=(2,2)).plot()
             ______________
            /              \
            \___  _  _*_  _/
                \/ \/   \/
             ___/\_/\___/\_
            /              \
            \_•____________/
        """
        
        if plat is None:
            super().__init__(crossinglist, name=name)
        else:
            super().__init__(plat.crossings(), name=plat.name())
        self.basepoint_bullet=basepoint_bullet
        self.basepoint_asterisk=basepoint_asterisk
        self.swap_basepoints=swap_basepoints
        c = len(self.__crossinglist__)
        nc = self.__ncusps__
        n = c + nc
        if basepoint_bullet[0] < 0 or basepoint_bullet[0] > c or basepoint_asterisk[0] < 0 or basepoint_asterisk[0] > c:
            raise ValueError('first coordinate of basepoint must be an integer in [0, num_crossings]')
        if basepoint_bullet[1] < 1 or basepoint_bullet[1] > 2*nc or basepoint_asterisk[1] < 1 or basepoint_asterisk[1] > 2*nc:
            raise ValueError('second coordinate of basepoint must be an integer in [1, 2*num_cusps]')

        lsft_varnames = ['t','u'] 
        lsft_varnames += ['q0'] if n==1 else ['q'+str(i) for i in range(n)]
        lsft_varnames += ['p0'] if n==1 else ['p'+str(i) for i in range(n)]
        self.__lsft_ring__ = PolynomialRing(GF(2), 2*n+2,lsft_varnames)

    def plot(self):
        """
        Print an ascii representation of the based plat.

        EXAMPLES::
        
            sage: K = LegendrianKnotTable['m4_1']; B = BasedPlat(K, basepoint_bullet=(3,2), basepoint_asterisk=(7,4)); B.plot()
             _____________________________
            /                             \
            \_  ____  _•_______________  _/
              \/    \/                 \/
             _/\____/\___  _______  ___/\_
            /            \/       \/      \
            \____  ______/\_  _  _/\_*____/
                 \/         \/ \/
             ____/\_________/\_/\_________
            /                             \
            \_____________________________/
        """

        self._build_based_plat_ascii()
        for line in self.__plat_ascii__:
            print(line)


    def hamiltonian(self):
        """
        Return the hamiltonian of this based plat.

        OUTPUT:
        
        An element of self.lsft_ring().

        EXAMPLES::
            sage: BasedPlat(LegendrianKnotTable['m3_1']).hamiltonian()
            t*q0*q1*q2*p3 + q0*q1*q2*p4 + t*q0*p3 + t*q2*p3 + p0*p1 + p1*p2 + q0*p4 + q2*p4 + p3 + p4
        """

        try:
            return self.__hamiltonian__
        except AttributeError:
            nc = self.ncusps()
            n = nc + self.ncrossings()
            alpha = [[0]*(2*nc+1) for i in range(0,2*nc+1)]
            hamiltonian = 0
            for i in range(nc):
                alpha[2*i+1][2*(i+1)] = 1
            bx = self.basepoint_bullet[0]
            bz = self.basepoint_bullet[1]
            t = self.lsft_ring()('t')
            u = self.lsft_ring()('u')
            if bx == 0:
                if bz % 2 == 0:
                    alpha[bz-1][bz] = u
                else:
                    alpha[bz][bz+1] = t
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

                if bx == a+1:
                    for i in range(1,bz):
                        alpha[i][bz] = alpha[i][bz]*u
                    for i in range(bz+1, 2*nc+1):
                        alpha[bz][i] = alpha[bz][i]*t
            
            #for row in alpha:
            #    print(row)
            #print('-------------')

            for i in range(nc):
                p = self.lsft_ring().gen(i+2+self.ncrossings()+n)
                hamiltonian = hamiltonian + p + p*alpha[2*i+1][2*i+2]
            
            self.__hamiltonian__ = hamiltonian
            return self.__hamiltonian__

    def lsft_ring(self):
        """
        Return the commutative lsft ring. 

        OUTPUT:

        The polynomial ring over Z_2 with variables p_i, q_i for each crossing and right cusp of the plat, as well as variables t, u.
        p_i corresponds to a positive end at the i^th crossing, while q_i corresponds to a negative end at the i^th crossing.

        EXAMPLES::

           sage: BasedPlat(crossinglist=[2,2,2]).lsft_ring()
           Multivariate Polynomial Ring in t, u, q0, q1, q2, q3, q4, p0, p1, p2, p3, p4 over Finite Field of size 2
        """
        return self.__lsft_ring__

    def lsft_differentials(self):
        """
        Return the lsft differentials of the generators of the lsft algebra associated to this based plat.

        OUTPUT:
        
        A dictionary whose keys are the generators of self.lsft_ring() and whose values are polynomials in self.lsft_ring().

        EXAMPLES::
            sage: BasedPlat(LegendrianKnotTable['m5_1']).lsft_differentials()
            {p0: t*q1*q2*q3*q4*p5 + q1*q2*q3*q4*p6 + t*q1*q2*p5 + t*q1*q4*p5 + t*q3*q4*p5 + q0*p0^2 + q1*p0*p1 + q2*p0*p2 + q3*p0*p3 + q4*p0*p4 + q1*q2*p6 + q1*q4*p6 + q3*q4*p6 + t*p5 + p6,
            q0: q0*q1*p1 + q0*q2*p2 + q0*q3*p3 + q0*q4*p4 + p1,
            p1: t*q0*q2*q3*q4*p5 + q0*q2*q3*q4*p6 + t*q0*q2*p5 + t*q0*q4*p5 + q0*p0*p1 + q1*p1^2 + q2*p1*p2 + q3*p1*p3 + q4*p1*p4 + q0*q2*p6 + q0*q4*p6,
            q1: q0*q1*p0 + q1*q2*p2 + q1*q3*p3 + q1*q4*p4 + p0 + p2,
            p2: t*q0*q1*q3*q4*p5 + q0*q1*q3*q4*p6 + t*q0*q1*p5 + t*q3*q4*p5 + q0*p0*p2 + q1*p1*p2 + q2*p2^2 + q3*p2*p3 + q4*p2*p4 + q0*q1*p6 + q3*q4*p6 + t*p5 + p6,
            q2: q0*q2*p0 + q1*q2*p1 + q2*q3*p3 + q2*q4*p4 + p1 + p3,
            p3: t*q0*q1*q2*q4*p5 + q0*q1*q2*q4*p6 + t*q0*q4*p5 + t*q2*q4*p5 + q0*p0*p3 + q1*p1*p3 + q2*p2*p3 + q3*p3^2 + q4*p3*p4 + q0*q4*p6 + q2*q4*p6,
            q3: q0*q3*p0 + q1*q3*p1 + q2*q3*p2 + q3*q4*p4 + p2 + p4,
            p4: t*q0*q1*q2*q3*p5 + q0*q1*q2*q3*p6 + t*q0*q1*p5 + t*q0*q3*p5 + t*q2*q3*p5 + q0*p0*p4 + q1*p1*p4 + q2*p2*p4 + q3*p3*p4 + q4*p4^2 + q0*q1*p6 + q0*q3*p6 + q2*q3*p6 + t*p5 + p6,
            q4: q0*q4*p0 + q1*q4*p1 + q2*q4*p2 + q3*q4*p3 + p3,
            p5: 0,
            q5: t*q0*q1*q2*q3*q4 + t*q0*q1*q2 + t*q0*q1*q4 + t*q0*q3*q4 + t*q2*q3*q4 + q5^2*p5 + t*q0 + t*q2 + t*q4 + 1,
            p6: 0,
            q6: q0*q1*q2*q3*q4 + q0*q1*q2 + q0*q1*q4 + q0*q3*q4 + q2*q3*q4 + q6^2*p6 + q0 + q2 + q4 + 1,
            t: 0,
            u: 0}
        """

        try:
            return self.__lsft_differentials__
        except AttributeError:
            self.__lsft_differentials__ = dict()
            hamiltonian = self.hamiltonian()
            string_differentials = self.string_differentials()
            n = self.ncusps()+self.ncrossings()
            for i in range(n):
                p_i = self.lsft_ring().gen(i+2+n)
                q_i = self.lsft_ring().gen(i+2)
                self.__lsft_differentials__[p_i] = self.bracket(hamiltonian, p_i) + string_differentials[p_i]
                self.__lsft_differentials__[q_i] = self.bracket(hamiltonian, q_i) + string_differentials[q_i]
            t = self.lsft_ring()('t')
            u = self.lsft_ring()('u')
            self.__lsft_differentials__[t] = 0
            self.__lsft_differentials__[u] = 0
            return self.__lsft_differentials__


    def string_differentials(self):
        """
        Return the string differentials of the generators of the lsft algebra associated to this based plat.

        OUTPUT:
        
        A dictionary whose keys are the generators of self.lsft_ring() and whose values are polynomials in self.lsft_ring().

        EXAMPLES::
            sage: BasedPlat(LegendrianKnotTable['m5_1']).string_differentials()
            {q0: q0*q1*p1 + q0*q2*p2 + q0*q3*p3 + q0*q4*p4,
            p0: q0*p0^2 + q1*p0*p1 + q2*p0*p2 + q3*p0*p3 + q4*p0*p4,
            q1: q0*q1*p0 + q1*q2*p2 + q1*q3*p3 + q1*q4*p4,
            p1: q0*p0*p1 + q1*p1^2 + q2*p1*p2 + q3*p1*p3 + q4*p1*p4,
            q2: q0*q2*p0 + q1*q2*p1 + q2*q3*p3 + q2*q4*p4,
            p2: q0*p0*p2 + q1*p1*p2 + q2*p2^2 + q3*p2*p3 + q4*p2*p4,
            q3: q0*q3*p0 + q1*q3*p1 + q2*q3*p2 + q3*q4*p4,
            p3: q0*p0*p3 + q1*p1*p3 + q2*p2*p3 + q3*p3^2 + q4*p3*p4,
            q4: q0*q4*p0 + q1*q4*p1 + q2*q4*p2 + q3*q4*p3,
            p4: q0*p0*p4 + q1*p1*p4 + q2*p2*p4 + q3*p3*p4 + q4*p4^2,
            q5: q5^2*p5,
            p5: 0,
            q6: q6^2*p6,
            p6: 0}
        """

        try:
            return self.__string_differentials__
        except AttributeError:
            n = self.ncrossings() + self.ncusps()

            betas = dict()
            holomorphic_corners = dict()

            for i in range(n):
                q = self.lsft_ring().gen(2+i)
                betas[q] = 0
                holomorphic_corners[q] = False

            start_x = self.basepoint_bullet[0]
            start_z = self.basepoint_bullet[1]

            #first go right
            self._traverse(betas, holomorphic_corners, start_x, start_z, True)
    
            #then go left
            self._traverse(betas, holomorphic_corners, start_x, start_z, False)

            self.__string_differentials__ = dict()
            for i in range(n):
                q = self.lsft_ring().gen(2+i)
                p = self.lsft_ring().gen(2+n+i)
                self.__string_differentials__[q] = betas[q]*q if holomorphic_corners[q] else betas[q]*q+q*p*q
                self.__string_differentials__[p] = betas[q]*p+p*q*p if holomorphic_corners[q] else betas[q]*p

            return self.__string_differentials__

    def bracket(self, x, y):
        """
        Return the SFT bracket of x and y.

        INPUT:

        - ``x`` -- an element of self.lsft_ring()

        - ``y`` -- an element of self.lsft_ring()

        OUTPUT:
        
        An element of self.lsft_ring().

        EXAMPLES::

            sage: B=BasedPlat(LegendrianKnotTable['m5_1'])
            sage: p0=B.lsft_ring()('p0')
            sage: q0=B.lsft_ring()('q0')
            sage: B.bracket(p0,q0)
            1
        """

        x_grad = x.gradient()
        y_grad = self.dual(y).gradient()
        bracket_xy = 0
        for i in range(len(x_grad)):
            bracket_xy += x_grad[i] * self.dual(y_grad[i])
        return bracket_xy


    def d(self, x):
        """
        Return the LSFT differential of an element of the algebra.

        INPUT:

        - ``x`` -- an element of self.lsft_ring()

        OUTPUT:
        
        An element of self.lsft_ring().

        EXAMPLES::

            sage: B=BasedPlat(LegendrianKnotTable['m3_1'])
            sage: x=B.lsft_ring()('p0*q0')
            sage: B.d(x)
            t*q0*q1*q2*p3 + q0^2*p0^2 + q0*q1*q2*p4 + t*q0*p3 + p0*p1 + q0*p4
        """

        dx = 0
        for i in range(0, self.lsft_ring().ngens()):
            gen = self.lsft_ring().gen(i)
            dx = dx + x.derivative(gen) * self.lsft_differentials()[gen]
        return dx

    def dual(self, x):
        """
        Return the dual of an element of the algebra.

        INPUT:

        - ``x`` -- an element of self.lsft_ring()

        OUTPUT:
        
        An element of self.lsft_ring().

        EXAMPLES::

            sage: B=BasedPlat(LegendrianKnotTable['m3_1'])
            sage: x=B.lsft_ring()('p0*q1')
            sage: B.dual(x)
            q0*p1
        """
        x_str = str(x)
        x_str = x_str.replace('p', 'w')
        x_str = x_str.replace('q', 'p')
        x_str = x_str.replace('w', 'q')
        return self.lsft_ring()(x_str)
    

    def _build_based_plat_ascii(self):
        try:
            return self.__plat_ascii__
        except AttributeError:
            self._build_left_cusps()
            self._build_basepoints(0)
            for i in range(self.ncrossings()):
                self._build_extension()
                self._build_crossing(self.__crossinglist__[i])
                self._build_basepoints(i+1)

            #self._build_basepoints(self.ncrossings()+1)
            self._build_extension()
            self._build_right_cusps()
        return self.__plat_ascii__

    def _build_basepoints(self, i):
        if self.basepoint_bullet[0] == i and self.basepoint_asterisk[0] == i and self.basepoint_bullet[1] == self.basepoint_asterisk[1]:
            self._build_extension()
            self._build_basepoint(self.basepoint_bullet[1], self.swap_basepoints)
            self._build_extension()
            self._build_basepoint(self.basepoint_bullet[1], not self.swap_basepoints)
        else:
            if self.basepoint_bullet[0] == i:
                self._build_extension()
                self._build_basepoint(self.basepoint_bullet[1], True)
            if self.basepoint_asterisk[0] == i:
                self._build_extension()
                self._build_basepoint(self.basepoint_asterisk[1], False)

    def _build_basepoint(self, i, bullet=True):
        for j in range(len(self.__plat_ascii__)):
            strand_char = '_' if j%2==0 else ' '
            if j==2*i-2:
                strand_char = '•' if bullet else '*'
            self.__plat_ascii__[j] = self.__plat_ascii__[j] + strand_char

    

    def _traverse(self, betas, holomorphic_corners, x, z, right=True):
        target_x = self.basepoint_asterisk[0]
        target_z = self.basepoint_asterisk[1]
        target_reached = (x == target_x) and (z == target_z) and (self.swap_basepoints == right)
        beta = 0
        while not target_reached:
            (beta, x, z, right, target_reached) = self._step(betas, holomorphic_corners, beta, x, z, target_x, target_z, right)
            #print('beta = {}, x = {}, z = {}, right = {}, target_reached = {}'.format(beta, x, z, right, target_reached))

    #return tuple (new_beta, new_x, new_z, new_right, target_reached)
    def _step(self, betas, holomorphic_corners, beta, x, z, target_x, target_z, right=True):
        #right cusp reached - update cusp and turn back.
        if x == self.ncrossings() and right:
            q = self.lsft_ring().gen(2+self.ncrossings()+(z-1)//2)
            p = self.lsft_ring().gen(2+self.ncrossings()+self.ncusps()+self.ncrossings()+(z-1)//2)
            betas[q] = q*p
            holomorphic_corners[q] = True
            new_z = self._new_z_cusp(z)
            return (beta, x, new_z, False, False)
        elif x == 0 and not right:
            new_z = self._new_z_cusp(z)
            target_reached = (x==target_x) and ((z==target_z) or (new_z==target_z))
            return (beta, x, new_z, True, target_reached)
        else:
            if right:
                crossing = self.__crossinglist__[x]
                new_x = x+1
                new_z = self._new_z_crossing(z, crossing)
                target_reached = (new_x == target_x) and (new_z == target_z)
                if z!=new_z:
                    q = self.lsft_ring().gen(2+x)
                    p = self.lsft_ring().gen(2+self.ncrossings()+self.ncusps()+x)
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
                        p = self.lsft_ring().gen(2+self.ncrossings()+self.ncusps()+new_x)
                        betas[q] = betas[q] + beta
                        beta = beta+q*p
                        holomorphic_corners[q] = not holomorphic_corners[q]
                    return (beta, new_x, new_z, False, False)

    def _new_z_crossing(self, z, crossing):
        return z+1 if z==crossing else (z-1 if z==crossing+1 else z)

    def _new_z_cusp(self, z):
        return (z-1) if z%2==0 else (z+1)
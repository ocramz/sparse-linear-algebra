module Data.Sparse.TestData where

import Data.Sparse.Common
import Data.Complex
import Data.VectorSpace
{-

example 0 : 2x2 linear system

[1 2] [2] = [8]
[3 4] [3]   [18]

[1 3] [2] = [11]
[2 4] [3]   [16]


-}


aa0 :: SpMatrix Double
aa0 = fromListDenseSM 2 [1,3,2,4]

-- b0, x0 : r.h.s and initial solution resp.
b0, x0, x0true, aa0tx0 :: SpVector Double
b0 = mkSpVR 2 [8,18]
x0 = mkSpVR 2 [0.3,1.4]


-- x0true : true solution
x0true = mkSpVR 2 [2,3]

aa0tx0 = mkSpVR 2 [11,16]







{- 4x4 system -}

aa1 :: SpMatrix Double
aa1 = sparsifySM $ fromListDenseSM 4 [1,0,0,0,2,5,0,10,3,6,8,11,4,7,9,12]

x1, b1 :: SpVector Double
x1 = mkSpVR 4 [1,2,3,4]

b1 = mkSpVR 4 [30,56,60,101]



{- 3x3 system -}
aa2 :: SpMatrix Double
aa2 = sparsifySM $ fromListDenseSM 3 [2, -1, 0, -1, 2, -1, 0, -1, 2]
x2, b2 :: SpVector Double
x2 = mkSpVR 3 [3,2,3]

b2 = mkSpVR 3 [4,-2,4]


aa22 = fromListDenseSM 2 [2,1,1,2] :: SpMatrix Double





{- 2x2 Complex system -}

aa0c :: SpMatrix (Complex Double)
aa0c = fromListDenseSM 2 [ 3 :+ 1, (-3) :+ 2, (-2) :+ (-1), 1 :+ (-2)]

b0c = mkSpVC 2 [3 :+ (-4), (-1) :+ 0.5]

x1c = mkSpVC 2 [2 :+ 2, 2 :+ 3]
b1c = mkSpVC 2 [4 :+ (-2), (-10) :+ 1]

aa2c :: SpMatrix (Complex Double)
aa2c = fromListDenseSM 2 [3, -3, -2, 1]








-- matlab : aa = [1, 2-j; 2+j, 1-j]
aa3c, aa3cx :: SpMatrix (Complex Double)
aa3c = fromListDenseSM 2 [1, 2 :+ 1, 2 :+ (-1), 1 :+ (-1)]

-- matlab : aaxaa = aa * aa
aa3cx = fromListDenseSM 2 [6, 5, 3 :+ (-4), 5:+ (-2)]







{-
matMat

[1, 2] [5, 6] = [19, 22]
[3, 4] [7, 8]   [43, 50]
-}

m1, m2, m1m2, m1', m2', m1m2', m2m1' :: SpMatrix Double
m1 = fromListDenseSM 2 [1,3,2,4]
m2 = fromListDenseSM 2 [5, 7, 6, 8]     
m1m2 = fromListDenseSM 2 [19, 43, 22, 50]

m1' = fromListSM (2,3) [(0,0,2), (1,0,3), (1,2,4), (1,2,1)]
m2' = fromListSM (3,2) [(0,0,5), (0,1,3), (2,1,4)]
m1m2' = fromListDenseSM 2 [10,15,6,13] 
m2m1' = fromListSM (3,3) [(0,0,19),(2,0,12),(0,2,3),(2,2,4)]

-- transposeSM
m1t :: SpMatrix Double
m1t = fromListDenseSM 2 [1,2,3,4]


--

{-
countSubdiagonalNZ
-}
m3 :: SpMatrix Double
m3 = fromListSM (3,3) [(0,2,3),(2,0,4),(1,1,3)] 






{- eigenvalues -}
aa3 :: SpMatrix Double
aa3 = fromListDenseSM 3 [1,1,3,2,2,2,3,1,1]

b3 = mkSpVR 3 [1,1,1] :: SpVector Double



-- aa4 : eigenvalues 1 (mult.=2) and -1
aa4 :: SpMatrix Double
aa4 = fromListDenseSM 3 [3,2,-2,2,2,-1,6,5,-4] 

aa4c :: SpMatrix (Complex Double)
aa4c = toC <$> aa4

b4 = fromListDenseSV 3 [-3,-3,-3] :: SpVector Double




aa5c :: SpMatrix (Complex Double)
aa5c = fromListDenseSM 4 cv where
  cv = zipWith (:+) [1..16] [16,15..1]






tm0, tm1, tm2, tm3, tm4, tm5, tm6 :: SpMatrix Double
tm0 = fromListSM (2,2) [(0,0,pi), (1,0,sqrt 2), (0,1, exp 1), (1,1,sqrt 5)]

tv0, tv1 :: SpVector Double
tv0 = mkSpVR 2 [5, 6]

tv1 = fromListSV 2 [(0,1)] 

-- wikipedia test matrix for Givens rotation

tm1 = sparsifySM $ fromListDenseSM 3 [6,5,0,5,1,4,0,4,3]

-- wp test matrix for QR decomposition via Givens rotation

tm2 = fromListDenseSM 3 [12, 6, -4, -51, 167, 24, 4, -68, -41]

-- λ> (l,u) <- lu tm2
-- λ> prd l

-- 1.00   , _      , _      
-- 0.50   , 1.00   , _      
-- -0.33  , 0.04   , 1.00   

-- λ> prd u

-- 12.00  , -51.00 , 4.00   
-- _      , 1.92e2 , -70.00 
-- _      , _      , -37.12 




tm3 = transposeSM $ fromListDenseSM 3 [1 .. 9]



--

tm4 = sparsifySM $ fromListDenseSM 4 [1,0,0,0,2,5,0,10,3,6,8,11,4,7,9,12]


tm5 = fromListDenseSM 3 [2, -4, -4, -1, 6, -2, -2, 3, 8] 


tm6 = fromListDenseSM 4 [1,3,4,2,2,5,2,10,3,6,8,11,4,7,9,12] 

tm7 :: SpMatrix Double
tm7 = a ^+^ b ^+^ c where
  n = 5
  a = mkSubDiagonal n 1 $ replicate n (-1)
  b = mkSubDiagonal n 0 $ replicate n 2
  c = mkSubDiagonal n (-1) $ replicate n (-1)

tvx7 = mkSpVR 5 [3,8,-12,4,9]

tvb7 = tm7 #> tvx7




tm8 :: SpMatrix Double
tm8 = fromListSM (2,2) [(0,0,1), (0,1,1), (1,1,1)]

tm8' :: SpMatrix Double
tm8' = fromListSM (2,2) [(0,0,1), (1,0,1), (1,1,1)]



tm9 :: SpMatrix Double
tm9 = fromListSM (4, 3) [(0,0,pi), (1,1, 3), (2,2,4), (3,2, 1), (3,1, 5)]





-- tvc0 <.> tvc1 = 5 
tvc0, tvc1, tvc2, tvc3 :: SpVector (Complex Double)
tvc0 = fromListSV 2 [(0,0), (1,2 :+ 1)]
tvc1 = fromListSV 2 [(0,0), (1, 2 :+ (-1))] 


-- dot([1+i, 2-i], [3-2i, 1+i]) = 2 - 2i
tvc2 = fromListDenseSV 2 [1 :+ 1,  2 :+ (-1)]
tvc3 = fromListDenseSV 2 [3 :+ (-2), 1 :+ 1 ]




-- Complex linear system

tmc4,tmc5, tmc6 :: SpMatrix (Complex Double)
-- tmc4: condition number = 4.4233
tmc4 = fromListDenseSM 3 [3:+1, 4:+(-1), (-5):+3, 2:+2, 3:+(-2), 5:+0.2, 7:+(-2), 9:+(-1), 2:+3]

-- tvc4 : unknown to be found
tvc4 = vc [1:+3,2:+2,1:+9]
-- bc4 : right-hand side
bc4 = tmc4 #> tvc4

tmc5 = fromListDenseSM 4 $ zipWith (:+) [16..31] [15,14..0]

tmc6 = fromListDenseSM 4 (zipWith (:+) [0..15] (replicate 16 1))





-- Rectangular real system

aa10 :: SpMatrix Double
aa10 = fromListDenseSM 3 [1,2,3,4,5,6]

x10, b10 :: SpVector Double
x10 = fromListDenseSV 2 [2,3]
b10 = aa10 #> x10

--

-- | Example 5.4.2 from G & VL
-- aa1 :: SpMatrix Double
-- aa1 = transpose $ fromListDenseSM 3 [1..12]

-- aa1 :: SpMatrix Double
-- aa1 = sparsifySM $ fromListDenseSM 4 [1,0,0,0,2,5,0,10,3,6,8,11,4,7,9,12]






-- l0 = [1,2,4,5,8]
-- l1 = [2,3,6]
-- l2 = [7]

-- v0,v1 :: V.Vector Int
-- v0 = V.fromList [0,1,2,5,6]
-- v1 = V.fromList [0,3,4,6]

-- -- e1, e2 :: V.Vector (Int, Double)
-- -- e1 = V.indexed $ V.fromList [1,0,0]
-- -- e2 = V.indexed $ V.fromList [0,1,0]

-- e1, e2:: CsrVector Double
-- e1 = fromListCV 4 [(0, 1)] 
-- e2 = fromListCV 4 [(1, 1)]
-- e3 = fromListCV 4 [(0, 1 :+ 2)] :: CsrVector (Complex Double)

-- e1c = V.indexed $ V.fromList [1,0,0] :: V.Vector (Int, Complex Double)

-- m0,m1,m2,m3 :: CsrMatrix Double
-- m0 = toCSR 2 2 $ V.fromList [(0,0, pi), (1,0,3), (1,1,2)]
-- m1 = toCSR 4 4 $ V.fromList [(0,0,1), (0,2,5), (1,0,2), (1,1,3), (2,0,4), (2,3,1), (3,2,2)]
-- m2 = toCSR 4 4 $ V.fromList [(0,0,1), (0,2,5), (2,0,4), (2,3,1), (3,2,2)]
-- m3 = toCSR 4 4 $ V.fromList [(1,0,5), (1,1,8), (2,2,3), (3,1,6)]

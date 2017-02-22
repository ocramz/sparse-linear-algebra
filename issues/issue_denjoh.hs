import Numeric.LinearAlgebra.Sparse
import Data.Sparse.Common

eMod = 210000.0


k :: Int -> Int -> Double -> Double -> Double
k 0 0 i l = 12* eMod * i /(l**3)
k 0 1 i l = 6 * eMod * i /(l**2)
k 0 2 i l = (- (k 0 0 i l))
k 0 3 i l = k 0 1 i l 
k 1 0 i l = k 0 1 i l 
k 1 1 i l = 4 * eMod *i / l
k 1 2 i l = - k 0 1 i l 
k 1 3 i l = 2* eMod * i /l
k 2 0 i l = - k 0 0 i l
k 2 1 i l = - k 0 1 i l
k 2 2 i l = k 0 0 i l
k 2 3 i l = - k 0 1 i l
k 3 0 i l = k 0 1 i l 
k 3 1 i l = k 1 3 i l
k 3 2 i l = - k 0 1 i l
k 3 3 i l = k 1 1 i l
k _ _ _ _ = error "Index out of bounds: Local Stiffness Matrix is only 4x4" 


localStiffness :: 
     [Double] -> [Double] -> Int -> [(Int, Int, Double)]
localStiffness [i1] [l1] 0 = [ (a , b, (k a b i1 l1) ) | a <- [0..3], b <- [0..3], not(a>1 && b> 1)]
localStiffness [i1,i2] [l1,l2] rn = 
    firstQ ++ rest
    where 
        firstQ = [ (rn + a, rn + b, ((k (a+2) (b+2) i1 l1) + (k a b i2 l2) )) | a <- [0..1], b <- [0..1]]
        rest = [ (rn + a, rn + b, ( k a b i2 l2 )) | a <- [0..3], b <- [0..3], not(a<2 && b<2)]
localStiffness (i1:i2:it) (l1:l2:lt) rn = 
    firstQ ++ sndTrdQ
    where 
        firstQ = [ (rn + a, rn + b, ((k (a+2) (b+2) i1 l1) + (k a b i2 l2) )) | a <- [0..1], b <- [0..1]]
        sndTrdQ = [ (rn + a, rn + b, ( k a b i2 l2 )) | a <- [0..3], b <- [0..3], not(a<2 && b<2), not(b>1 && a >1)]


createStiffnessMatrix :: [Double]
     -> [Double] -> [(Int, Int, Double)] -> Int -> SpMatrix Double
createStiffnessMatrix (i1:i2:iTail) (l1:l2:lTail) [] _ = 
   createStiffnessMatrix (i1:i2:iTail) (l1:l2:lTail) (localStiffness [i1] [l1] 0) 2
createStiffnessMatrix [i1, i2] [l1, l2] matrix rn =
    fromListSM (rn+4,rn+4) ( matrix ++ localStiffness [i1,i2] [l1,l2] (rn))
createStiffnessMatrix (i1:iTail) (l1:lTail) matrix rn = 
   createStiffnessMatrix iTail lTail ((localStiffness (i1:iTail) (l1:lTail) rn)++ matrix) (rn+2)

   
insertSprings :: (Foldable t, Num a) =>
     SpMatrix a -> t (IxRow, IxCol, a) -> SpMatrix a
insertSprings kMat springs = 
        foldl (addSpring) kMat springs
        where
                addSpring k (a, b, v) = insertSpMatrix a b ( ( k @@! (a,b) ) + v) k 


main = do
        let l=replicate 500 100
        let i=replicate 500 400000000
        let s = [(x,x,2000)|x<- [200,400..1000]]
        let k=createStiffnessMatrix i l [] 0
        let knew = insertSprings  k s
        let f= [ (x , 0) | x<-[0..9]] ++ [(1000,5000.0),(1001,819000000.0)] :: [(Int, Double)]
        let fVec = (fromListSV 1002 f) 
        --let res= (\(lo, up)-> luSolve lo up (dropSV 2 fVec)) $ lu (extractSubmatrixSM (\x -> x-2) (\x -> x-2) knew (2,1001) (2,1001))
        res <- extractSubmatrixSM (subtract 2) (subtract 2) knew (2,1001) (2,1001) <\> dropSV 2 fVec
        let myres = map (`lookupSV` res) [198,398..998]
        print myres

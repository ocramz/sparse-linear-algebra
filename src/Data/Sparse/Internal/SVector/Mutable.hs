{-# language FlexibleContexts, TypeFamilies #-}
module Data.Sparse.Internal.SVector.Mutable where

import qualified Data.Vector as V 
import qualified Data.Vector.Mutable as VM

import Control.Arrow ((&&&))
import Control.Monad.Primitive

data SMVector m a = SMV { smvDim :: {-# UNPACK #-} !Int,
                          smvIx :: V.Vector Int,
                          smvVal :: VM.MVector (PrimState m) a }


-- instance Show a => Show (SMVector m a) where
--   show (SMV n ix v) = unwords ["SMV (",show n,"),",show nz,"NZ:",show (V.zip ix v)]
--     where nz = V.length ix

-- instance Functor (SMVector m) where
--   fmap f (SMV n ix v) = SMV n ix (fmap f v)

-- instance Foldable (SMVector m) where
--   foldr f z (SMV _ _ v) = foldr f z v

-- instance Traversable (SMVector m) where
--   traverse f (SMV n ix v) = SMV n ix <$> traverse f v



-- | traverse the sparse mutable vector and operate on the union set between the indices of the immutable vector `v` and those of the mutable one `vm`, _overwriting_ the values in vm.
-- Invariant: the index set of `v` is a strict subset of that of `vm` (i.e. we assume that `vm` is preallocated properly, or, we assume there won't be any out-of-bound writes attempted)
unionWithM_prealloc g z v vm@(SMV nvm vmix vmv) = undefined








-- | unionWithM takes the union of two sparse mutable vectors given a binary function and a neutral element.
unionWithM :: PrimMonad m =>
     (a -> a -> b)
     -> a
     -> SMVector m a
     -> SMVector m a
     -> m (VM.MVector (PrimState m) (Int, b))
unionWithM g z (SMV n1 ixu uu) (SMV n2 ixv vv) = do
  vm <- VM.new n
  (vm', nfin) <- go ixu uu ixv vv 0 vm
  let vm'' = VM.take nfin vm'
  return vm''
  where
    headTail = V.head &&& V.tail
    n = min n1 n2
    go iu u_ iv v_ i vm
          | (VM.null u_ && VM.null v_) || i == n = return (vm , i)
          | VM.null u_ = do
              v0 <- VM.read v_ 0
              VM.write vm i (V.head iv, g z v0)
              go iu u_ (V.tail iv) (VM.tail v_) (i + 1) vm
          | VM.null v_ = do
              u0 <- VM.read u_ 0
              VM.write vm i (V.head iu, g u0 z)
              go (V.tail iu) (VM.tail u_) iv v_ (i + 1) vm
          | otherwise =  do
             u <- VM.read u_ 0
             v <- VM.read v_ 0
             let us = VM.tail u_
                 vs = VM.tail v_
             let (iu1, ius) = headTail iu
                 (iv1, ivs) = headTail iv
             if iu1 == iv1 then do VM.write vm i (iu1, g u v)
                                   go ius us ivs vs (i + 1) vm
                           else if iu1 < iv1 then do VM.write vm i (iu1, g u z)
                                                     go ius us iv v_ (i + 1) vm
                                             else do VM.write vm i (iv1, g z v)
                                                     go iu u_ ivs vs (i + 1) vm


-- -- * helpers

-- -- both :: Arrow arr => arr b c -> arr (b, b) (c, c)
-- -- both f = f *** f



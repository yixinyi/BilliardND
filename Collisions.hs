module Collisions
( newPoints
, compareWith
)
where 

import Geometry
import Data.List (takeWhile, dropWhile)
import Control.Monad (join)


{- 
   The trajectory of a point-like particle within a limited region (i.e. a billiard)
   is characterized by an ordered set of N points = {p_1, p_2, ..., p_N}, that we call collisions,
   where each point has a coordinate p_i = (p1_i, ..., pn_i), n = dimension of space
   Let
        distPoints = {d(1,2), d(2,3), ..., d(n-1,n)}
        cumDistPoints = {d(1,2), d(1,2) + d(2,3), d(1,2) + d(2,3) + d(3,4), ... }
   where d(i,j) is the distance between points p_i and p_j
-}

distPoints :: [Coordinate] -> [Double]
distPoints ps = map2 [] distance ps

-- map2 applies a function f(x, y) to a list
map2 acc f [] = reverse acc
map2 acc f (x:xs) = if xs == []
                      then reverse acc
                      else map2 (new:acc) f xs
                    where new = f x (head xs)

cumDistPoints :: [Coordinate] -> [Double]
cumDistPoints ps = tail$scanl (+) 0 (distPoints ps) 

--undoCum [] = []
--undoCum l = (head l):(map2 [] (-) (reverse l))
                            
                            
                              
---------------------------------------------------------
-- Compare two different series --
---------------------------------------------------------
                    
{- 
   Given collisions: p = {p_0, p_1, p_2, ..., p_N}
                     q = {1_0, q_1, q_2, ..., q_N}
   Assume the same velocity for particle p and q
   
   We want to compare 2 trajectories at the same time.
   
   At the time when p_i happened, where was the particle q in its trajectory?
         p_i has travelled lp_i distance respect p_0
         If lq_{k-1} < lp_i < lq_k for some k, then return delta_{ki} = lq_k - lp_i
         There could be many lp_i satisfying the inequality above for one k.
     
-}

deltas :: [Double] -> [Double] -> [[Double]]
deltas lp lq 
    = lessThan [] lp lq
    where lessThan acc x [] = reverse acc -- points of x larger than y are ignored
          lessThan acc x (y:ys) = lessThan (dless:acc) rest ys
                                  where less = takeWhile ( <= y) x
                                        rest = dropWhile ( <= y) x
                                        dless = map (y-) less


-- Return a lists of new points q' on the trajectory of q at the same time when particle p collided with the walls
newPoints :: [Coordinate] -> [Coordinate] -> [Coordinate]
newPoints p q 
  = (head q):(join (helper [] q d us))
  where lp = cumDistPoints p
        lq = cumDistPoints q 
        -- unit vectors from the neighbouring points of the trajectory q
        -- {u_{0,1}, u_{1,2}, ...}
        us = map2 [] unitV q
        -- [[delta_{11}, delta_{12}],[delta_{21}],...] -- List of lists
        d = deltas lp lq
        -- q_i' = q_{k} - u_{k-1, k} * delta_{ki}
        helper acc qs ds [] = reverse acc
        helper acc (q:qs) (d:ds) (u:us) = helper (qnew:acc) qs ds us
                                        where new = [map (*x) y | x <- d, y <- [u]]
                                              qnew = map (zipWith (-) (head qs)) new
                                              
                                              
                                              
-- compare 2 collision sets with a function f
-- examples of f: 
--      \a b -> abs (a-b)
--      distance
compareWith f ps qs 
    = zip lp devs
    where qs' = newPoints ps qs
          devs = zipWith f ps qs'
          lp = cumDistPoints ps    -- acts as time

      


module SinaiND
( boxBilliard
, sinaiBilliard
, sinaiEBilliard
, history
)
where 

import Geometry
import Data.List (takeWhile, dropWhile)
import Control.Monad (join)

--------------------------------------------------
--------------------------------------------------
{- 
   Hypercube of dimension n
   The number of walls is 2*n
   The number of vertices is 2^n
   See more: https://en.wikipedia.org/wiki/Hypercube

   n-rectangle centered at origin, where boundaries are fixed by size = [bx, by, ...]
   Example:
         (-bx, by)_________(bx, by)  
              |              | 
              |              |
              |______________|
         (-bx, -by)        (bx, -by)
-}

boxBilliard size Line {point = p, direction = d}
    | nullQ d   = Line {point = p, direction = d} -- not moving
    | otherwise = Line {point = pSol, direction = dSol}
                where {- Function to label walls in the way below:
                            wall_i+ = Plane {pointPlane = [0,..., b_i, ..., 0], normal = [0,...,1,...,0]}
                            wall_i- = Plane {pointPlane = [0,..., -b_i, ..., 0], normal = [0,...,-1,...,0]}
                      -}                      
                      toWall i sign boundary = Plane {pointPlane = x, normal = n}
                                             where x = [if x == i then sign * boundary else 0 | x <- [1 .. dim]]  
                                                   n = [if x == i then sign else 0 | x <- [1 .. dim]] 
                                                   dim = length size
                                                               
                      -- Walls that the trajectory might intercept
                      walls = loop [] toWall 1 signs size
                            where loop acc f i s [] = reverse acc 
                                  loop acc f i (s:ss) (b:bs) = loop ((f i s b):acc) f (i+1) ss bs 
                                  signs = signsV d -- orientation
                                                               
                      -- Intersections of the trajectory with the walls -- a list of Maybe type
                      intersections = map (intersectLinePlane Line {point = p, direction = d}) walls 
                      
                      -- Distances of the intersections with respect to the origin -- a list of Maybe type
                      dists = map (fmap moduleV) intersections 
                      
                      -- The minimal distance, in order to determine the closest wall(s) the particle collides
                      minDist = minimum list
                              where list = fmap unJust [x | x <- dists, x /= Nothing] -- filter out Nothing
                              
                      -- The closest walls labelled by the positions, e.g. wall_1+ and wall_1- are both labeled by 1 
                      positionsWalls = positions [] (Just minDist) 1 dists
                      
                      -- Next location is the intersection with the closest wall(s) 
                      pSol = unJust$head$drop i intersections
                           where i = (head positionsWalls) - 1 
                           
                      -- Next direction
                      dSol = reflectionWall positionsWalls d 1

unJust (Just a) = a 

-- Positions of x in the list (d:ds)
positions acc x i [] = reverse acc
positions acc x i (d:ds) = if x == d
                           then positions (i:acc) x (i+1) ds
                           else positions acc x (i+1) ds
                           
---- Not tail-recursive version
--positions x i [] = [] 
--positions x i (d:ds) = if x == d
--                       then i:(positions x (i+1) ds)
--                       else positions x (i+1) ds

{- 
   Next direction: multiply (-1) in directions when a wall was hit, 
   i.e. if wall_1+ was hit
        then the next direction is (-1, 0, ..., 0).direction
   i.e. if [wall_1+, wall_3-] were hit (so, an edge where the walls cross)
        then the next direction is (-1, 0, -1, ..., 0).direction
-}

-- direction = d:ds
reflectionWall _ [] i = []
reflectionWall [] ds i = ds
reflectionWall (p:ps) (d:ds) i = if i == p
                                then -d:(reflectionWall ps ds (i+1))
                                else d:reflectionWall (p:ps) ds (i+1)



---------------------------------------------------------------------------------
-- n-Sinai billiard = n-rectangle + unit n-sphere, both centered at the origin -- 

sinaiBilliard size l
    = if 1 + e > onSphere -- due to numerical precision, we cannot use equality
      then boxBilliard size l
      else case intersectLineSphere l sphere of 
              Nothing -> boxBilliard size l
              Just (x:xs) -> if xs == []
                             then boxBilliard size l
                             else if distance2 p x < distance2 p (head xs)
                                 then Line {point = x, direction = reflectionSphere x d}
                                 else Line {point = head xs, direction = reflectionSphere (head xs) d}
    where onSphere = moduleV p
          e = 1e-8
          p = point l
          d = direction l
          dim = length size
          sphere = Sphere {origin = [0 | x <- [1..dim]] , radius = 1} 
          
          

-- n-SinaiE billiard = n-rectangle + n-ellipsoid, both centered at the origin 

sinaiEBilliard size ellipsoid l
    = case intersectLineEllipsoid l ellipsoid of 
          Nothing -> boxBilliard size l
          Just (x:xs) -> if xs == [] 
                         then boxBilliard size l
                         else if e > px -- on the ellipsoid
                              then boxBilliard size l
                              else if px < distance2 p (head xs)
                                   then Line {point = x, direction = reflectionEllipsoid ellipsoid (normalTPlane x) d}
                                   else Line {point = head xs, direction = reflectionEllipsoid ellipsoid (normalTPlane (head xs)) d}
                         where e = 1e-8 
                               px = distance2 p x   
                               p = point l
                               d = direction l
                               normalTPlane p = zipWith (*) p coefs -- normal of the tangent plane   
                               coefs = coefficients ellipsoid   

{- 
   Compute the reflected direction:
     1. find the plane that is expanded by the incident vector i and the normal of the tangent plane p
     2. find a orthonormal basis on this plane
        one basis is the normalized incident vector i
        the other basis is determined by:
              q = s p + t i
              s and t are parameters that we fix by demanding:
                  q.p = 0 and q.q = 1
              we obtain:
                  t = - s * |p| / (cos alpha), alpha is the angle between p and i
                  s = +- cos alpha / (|p| * sin alpha )
              we choose the positive sign:    
                  q = + (cos alpha u_p - u_i) / (sin alpha)
              
              2 scenarios:
                          
                        / r
                       / 
              ===p====>---------- 
              |        \
              | q_+     \ i
              v

              ex: p = (0,1), i = (-1, 1)/ (sqrt 2), alpha = pi/2 + pi /4
                  cos alpha = - sin alpha = 1/(sqrt 2)
                  we get:
                  q_+- = -+ (1, 0)

              ^ q_+
              |         / i
              |        / 
              ===p====>---------- 
                       \
                        \ r
              
              ex: p = (0,1), i = (-1, -1)/ (sqrt 2), alpha = -(pi/2 + pi /4)
                  cos alpha = sin alpha = - 1/(sqrt 2)
                  we get:
                  q_+- = +- (1, 0)
                  
      3. the reflected vector in terms of the unit vectors u_p and u_q
         r = |i| (- cos alpha u_p - sin alpha u_q)
         
-}


reflectionSphere p i -- p lies on the unit sphere, hence |p| = 1
    = if parallelQ p i
      then map (* (-1)) i
      else r
    where r = zipWith (+) r1 r2
          r1 = map (* (- cosAngle * imodule)) p
          r2 = map (* (- sinAngle * imodule)) q
          q = zipWith (\u v -> (cosAngle * u - v / imodule) / sinAngle) p i
          imodule = moduleV i
          cosAngle = cosineVs p i
          sinAngle = sqrt (1 - cosAngle^2)
          

reflectionEllipsoid e p i -- p is the normal vector, i is the incident vector
    = if parallelQ p i
      then map (* (-1)) i
      else r
    where r = zipWith (+) r1 r2 -- reflected 
          r1 = map (* (- cosAngle * imodule)) up
          r2 = map (* (- sinAngle * imodule)) uq
          imodule = moduleV i
          up = map (/ (moduleV p)) p -- unit vector at direction to point p
          uq = zipWith (\u v -> (cosAngle * u - v / imodule) / sinAngle) up i           
          cosAngle = cosineVs p i
          sinAngle = sqrt (1 - cosAngle^2)
          
          
          
          

--------------------------------------------           
-- history prints out a list of locations --
--------------------------------------------

history n billiard line = map point collisions
                        where collisions = take n (iterate billiard line)                            


-- Testing
sinai2 = history 20 (sinaiBilliard [2, 2]) Line {point = [2, 0], direction = [-1,0.1]}
    





                
                                        

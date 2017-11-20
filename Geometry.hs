module Geometry
( Line(..)
, Ellipsoid(..)
, Sphere(..)
, Plane(..)
, Coordinate
, Vector
, signsV
, unitV
, distance
, distance2
, moduleV
, module2V
, scalarProduct
, cosineVs
, nullQ
, orthogonalQ
, parallelQ
, semiaxis
, intersectLinePlane
, intersectLineSphere
, intersectLineEllipsoid
) where
                                     
type Coordinate = [Double]
type Vector = [Double]


data Line = Line { point :: Coordinate
                 , direction :: Vector
                 } deriving (Eq, Show)
-- A line in any dimension is described parametrically:
-- r = r0 + d t
-- where t is the parameter and the rest are vectors.   


data Ellipsoid = Ellipsoid { originEllipsoid :: Coordinate
                           , coefficients :: [Double]
                           } deriving (Eq, Show)  
-- Equation of an ellipsoid:
-- c1 (x - x0)^2 + ... + c2 (z - z0)^2 = 1      
-- c_i = 1 / a_i^2, where a_i are the semiaxis 
                 
data Sphere = Sphere { origin :: Coordinate
                     , radius :: Double
                     } deriving (Eq, Show)  
-- Equation of a hypersphere (a special case of an ellipsoid):
-- (r - r0)^2 = 0 -> (x - x0)^2 + ... + (z - z0)^2 = r^2                     

data Plane = Plane { pointPlane :: Coordinate
                   , normal :: Vector
                    } deriving (Eq, Show) 
-- Equation of a plane is determined by the scalar product of the normal vector and the position vector minus a point on the plane:
-- n . (r - r0) = 0     


--------------------------------------------------                                    
-- semiaxis 
semiaxis:: Ellipsoid -> Coordinate
semiaxis e = map (\c -> 1/(sqrt c)) coef
                   where coef = coefficients e                                    
                                    
--------------------------------------------------
-- distance between two points
distance p q = sqrt (distance2 p q)

-- distance squared
distance2 p q = sum$zipWith (\a b -> (a-b)^2) p q


--------------------------------------------------
-- null vector test
nullQ v = if elem False bools
          then False
          else True
        where bools = map (== 0) v

-- sign of each component of a vector
signsV v =  [if x < 0 then -1 else 1 | x <- v]

-- unit vector from u to v
unitV u v = map (/(moduleV diff)) diff
          where diff = zipWith (-) v u
             
-- module of a vector
moduleV v = sqrt$module2V v

-- module squared of a vector 
module2V v = sum [x^2 | x <- v]

-- scalar product
scalarProduct u v = sum$zipWith (*) u v

-- cosine of the angle between 2 vectors
cosineVs u v = (scalarProduct u v) / moduleV u / moduleV v 

-- orthogonal test for a given numerical precision
orthogonalQ u v = if e > zero && zero > - e
                then True
                else False
              where zero = abs$cosineVs u v
                    e = 1e-8
                    
-- (anti-)parallel test for a given numerical precision
parallelQ u v = if 1 + e > one && one > 1 - e
                then True
                else False
              where one = abs (cosineVs u v)
                    e = 1e-8
                                                 
----------------------------------------------------
---- intersection between a line and a plane
intersectLinePlane Line {point = p, direction = d} Plane {pointPlane = x, normal = n} 
    | orthogonalQ n d = Nothing -- pararell lines (including lines on the plane) 
    | otherwise       = Just inters
                      where t = - (scalarProduct n px) / (scalarProduct n d)
                            px = zipWith (-) p x
                            inters = zipWith (+) p (map (*t) d) 
                              
---- intersection between a line and a sphere
---- Solve for a + b t + c t^2 = 0
intersectLineSphere Line {point = p, direction = d} Sphere {origin = o, radius = r}
    | discriminant <  0 = Nothing
    | discriminant == 0 = Just [inters t0]
    | otherwise         = Just [inters tp, inters tm]
                        where po = zipWith (-) p o
                              a = module2V po - r^2
                              b = 2 * (sum$zipWith (*) d po)
                              c = module2V d
                              discriminant = b^2 - 4 * a * c
                              t0 = -b / 2 / c
                              tp = t0 + sqrt discriminant / 2 / c
                              tm = t0 - sqrt discriminant / 2 / c
                              inters t = zipWith (+) po (map (*t) d)      
                                                      
---- intersection between a line and an ellipsoid
---- Solve for a + b t + c t^2 = 0
intersectLineEllipsoid Line {point = p, direction = d} Ellipsoid {originEllipsoid = o, coefficients = coef}
    | discriminant <  0 = Nothing
    | discriminant == 0 = Just [inters t0]
    | otherwise         = Just [inters tp, inters tm]
                        where po = zipWith (-) p o
                              a = (sum$zipWith (\x y -> x * y^2) coef po) - 1
                              b = 2 * (sum$zipWith (*) coef (zipWith (*) d po))
                              c = sum$zipWith (\x y -> x * y^2) coef d
                              discriminant = b^2 - 4 * a * c
                              t0 = -b / 2 / c
                              tp = t0 + sqrt discriminant / 2 / c
                              tm = t0 - sqrt discriminant / 2 / c
                              inters t = zipWith (+) po (map (*t) d)                                                        
                              

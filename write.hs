import System.IO
import SinaiND
import Geometry
import Collisions

main = do
     let n = 3
         boundaries = [1 | x <- [1..n]]
         e = 1e-4
         length = 30
         dim = round n -- to get integral
         
         p = history length (sinaiBilliard boundaries) Line {point = [1 | x <- [1..n]], direction = [-x/n | x <- [1..n]]}
         q = history length (sinaiBilliard boundaries) Line {point = [1 | x <- [1..n]], direction = [-x/n-e | x <- [1..n]]}
                
     writeFile ("results/p"++ (show dim) ++ ".txt") (show boundaries ++ "\n" ++ show p)
     writeFile ("results/q" ++ (show dim) ++ ".txt") (show boundaries ++ "\n" ++ show q)
     writeFile ("results/qNew" ++ (show dim) ++ ".txt") (show (newPoints p q))
     writeFile ("results/deviations" ++ (show dim) ++ ".txt") (show (compareWith distance p q))     










        
    

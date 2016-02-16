import Interface
import Snp_parser
import Snp_windows

import System.Environment
import Control.Monad (liftM)
import Control.Monad.Trans.Except
import Data.Either
import Data.List
import qualified Data.Text.Lazy as T
import qualified Data.Text.Lazy.IO as TIO

main :: IO ()
main = getOptions >>= print


type Chrom = T.Text

data WindowStats = WindowStats {
    samples :: Int
  , mean    :: Double
  , start   :: Int
  , end     :: Int
  , chrom   :: Chrom
}

instance Show WindowStats where
  show (WindowStats s m st e c) =
     show c ++ "\t" ++ show st ++ "\t" ++ show e ++ "\t" ++ show s ++"\t"++ show m

meanWindow :: Chrom -> Window Double -> WindowStats
meanWindow c w = WindowStats (round n) (sum / n) (Snp_windows.start w) (Snp_windows.end w) c
  where
    (n, sum) = foldl' (\(n', s) x -> (n' + 1, s + x)) (0.0,0.0) (wData w)

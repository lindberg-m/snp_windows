import Interface
import Snp_parser
import Calc_windows

import System.IO                   (hGetContents, stdin)
import System.Exit                 (die)
import Control.Monad               (liftM)
import Control.Monad.Reader        (runReader)
import Control.Monad.Trans.Except  (runExcept)
import Data.Monoid
import Data.Either                 (rights, isRight)
import Data.Text.Lazy              (unpack)
import qualified Data.Text.Lazy.IO as TIO

main :: IO ()
main = do 
  opts <- getOptions 
  let getContent = case input opts of
         Nothing  -> TIO.hGetContents stdin
         (Just x) -> TIO.readFile x
      parser = if nohead opts
               then parseSNPsFromLineNr 1
               else parseSNPsFromLineNr 2
      config = WindowConfig (wSize opts) (wStep opts)

  snps <- liftM parser getContent
  printResults  $ getSNPs config snps

-- (list of chromosomes-windows, failed parsings (starting with a left value))
type SNPparseResults a b = ([(Chrom, [Window b])], [Either String (SNP a)])

printResults :: SNPparseResults Double (Sum Double) -> IO ()
printResults (res, failed) = do
   mapM_ putStrLn . concat $ map getRes res
   mapM_ catchError failed
   where
     catchError (Left err) = die $ "\nError: \n " ++ err ++ "\n"
     catchError _          = return ()

     getRes (chr, w)       = map (showWindow chr) w
     showWindow chr wind   = unpack chr ++ "\t" ++ show (start wind) ++ "\t" ++ show (end wind) ++
                             "\t" ++ (show . windowSamples $ wData wind) ++ "\t" ++ showMean
       where showMean = show (sum / fromIntegral n)
             n        = windowSamples $ wData wind
             sum      = getSum . windowDat $ wData wind
                       

{- 
 - Provided a config for window-sizes, generate a tuple
 - with all windows in as first element, and snps
 - not mapped to windows in the second element - starting
 - with the first error code from the snp-parsing
 -}
getSNPs :: WindowConfig -> [ParsedSNP Double] -> SNPparseResults Double (Sum Double)
getSNPs wc = mapFirst ((flip runReader wc) . calcWindows) . span isRight . map runExcept
  where
    mapFirst f (x,y) = (f $ rights x, y)


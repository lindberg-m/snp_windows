module Snp_parser (parseSNPs, parseSNPs', ParseError (..), ParsedSNP) where

import Snp_windows

import Control.Applicative
import System.IO 
import System.Exit                   (die, exitFailure)
import Control.Monad.Trans           (lift, liftIO)
import Control.Monad.Trans.Except    (Except, ExceptT, ExceptT (..), runExceptT, catchE, throwE)
import Data.Text.Lazy                (Text, pack, unpack, splitOn)
import Data.Text.Lazy.Read           (decimal, double)

data ParseError = UnsortedSNPs 
                | PositionUnparsed String
                | ValueUnparsed String
                | TooFewFields 
                | ParseError String
                deriving (Show, Eq)

type LineNumber = Int

type ParsedSNP a = Except String (SNP a)

{- Simplest and primary way of parsing a list of SNPs
 - is the function "parseSNPs"
 - The error (Left) value contains information about
 - what caused the failure as well as line number 
 -
 - Remaining parsers are more general 
 - (parseSNPs' & parseSNPs'')
 -
 -}
parseSNPs :: [Text] -> [ParsedSNP Double]
parseSNPs = parseSNPs'

parseSNPs' :: Monad m => [Text] -> [ExceptT String m (SNP Double)]
parseSNPs' = map parseErrorLineNumber .
             zip [1..] . 
             parseSNPs'' getDouble
  where
    parseErrorLineNumber (i,x) = catchE x (\e -> throwE $ errorLineMsg i e)
    getDouble = (fmap fst) . double . head

parseSNPs'' :: Monad m => ([Text] -> Either String a) -> [Text] -> [ExceptT ParseError m (SNP a)]
parseSNPs'' f txt = assertSNPorder $ map (toSnp . sepFields) txt
  where sepFields = splitOn $ pack "\t"
        toSnp (x1:x2:xs) = catchE (ExceptT . return $ decimal x2) (\e -> throwE $ PositionUnparsed e) >>=
                           \a -> catchE (ExceptT . return $ f xs) (\e -> throwE $ ValueUnparsed e) >>=
                           \b -> return $ SNP x1 (fst a) b
        toSnp _          = throwE TooFewFields

assertSNPorder :: Monad m => [ExceptT ParseError m (SNP a)] -> [ExceptT ParseError m (SNP a)] 
assertSNPorder = (foldr f []) . (map runExceptT)
  where
  f x l@(y:_) = let
    v = g <$> x <*> (runExceptT y)
    in (ExceptT v) : l
  f x [] = [ExceptT x]
  
  g a b = case (a,b) of
    (Right x, Right y) -> if pos x > pos y
                          then Left UnsortedSNPs
                          else a
    _                  -> a
      
errorLineMsg :: LineNumber -> ParseError -> String
errorLineMsg i UnsortedSNPs = "The SNPs are not in order (line " ++ show i ++
                              " and " ++ show (i + 1) ++ ")."
errorLineMsg i (PositionUnparsed e) = "The position-field could not be parsed at line " ++ 
                                       show i ++ ": " ++ e
errorLineMsg i (ValueUnparsed e) = "Could not parse the value-field(s) at line " ++ show i ++": " ++ e
errorLineMsg i TooFewFields  = "Missing fields at line " ++ show i ++ "." 
errorLineMsg i (ParseError e) = "Couldn't parse line " ++ show i ++ ": " ++ e

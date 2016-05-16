module Snp_parser (parseSNPsFromLineNr, 
                   ParseError (..), 
                   ParsedSNP) where

import Calc_windows                  (SNP(..))
import Control.Monad.Trans.Except    (Except, ExceptT, ExceptT (..), runExceptT, catchE, throwE)
import qualified Data.Text.Lazy      as T
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
 -}

parseSNPsFromLineNr :: Monad m => LineNumber -> T.Text -> [ExceptT String m (SNP Double)]
parseSNPsFromLineNr i = map catchErrorLineNumber .
                         enumerateFrom i .
                         drop (i - 1) .
                         parseSNPs getDouble
                         
  where enumerateFrom n = zip [n..] 
        catchErrorLineNumber (n,x) = catchE x (\e -> throwE $ errorLineMsg n e)
        getDouble = (fmap fst) . double . head


parseSNPs :: Monad m => ([T.Text] -> Either String a) -> T.Text -> [ExceptT ParseError m (SNP a)]
parseSNPs f = assertSNPorder . map (toSnp . sepFields) . T.lines
  where sepFields = T.splitOn $ T.pack "\t"
        toSnp (_:_:[])   = throwE TooFewFields
        toSnp (x1:x2:xs) = catchE (ExceptT . return $ decimal x2) (\e -> throwE $ PositionUnparsed e) >>=
                           \(a,_) -> catchE (ExceptT . return $ f xs) (\e -> throwE $ ValueUnparsed e) >>=
                           \b     -> return $ SNP x1 a b
        toSnp _          = throwE TooFewFields

assertSNPorder :: Monad m => [ExceptT ParseError m (SNP a)] -> [ExceptT ParseError m (SNP a)] 
assertSNPorder = (foldr f []) . (map runExceptT)
  where f x l@(y:_) = let
          v = g <$> x <*> (runExceptT y)
          in (ExceptT v) : l
        f x [] = [ExceptT x]
    
        g a b = case (a,b) of
          (Right x, Right y) -> if pos x > pos y && chrom x == chrom y
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

def gen_limit(alps):
    return alps.generation >= alps.max_generations

best_score = float('inf')
def conv_test(alps, threshold=0.05):
    global best_score
    delta = abs(best_score - alps.best.score)
    best_score = min(best_score, alps.best.score)
    return delta < threshold*best_score

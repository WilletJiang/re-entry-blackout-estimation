# The Way of Code
## *An Aesthetic Manifesto*

---

## Prefer **Composition** Over **Inheritance**

Inheritance is a promise you cannot keep.  
It binds you to ancestors you did not choose,  
carries debts you did not incur,  
and breaks when the world shifts beneath it.

**Composition is freedom.**  
Small things, combined. Each thing does one thing.  
You build castles from pebbles, not monoliths from mountains.

```
# Inheritance: a chain of assumptions
class A inherits from B inherits from C inherits from...
# One change in C, and everything shatters.

# Composition: a constellation of choices
f(g(h(x)))
# Each function stands alone. Each can be replaced.
```

> *Objects are data with behavior bolted on.  
> Functions are behavior with data threaded through.  
> Choose wisely.*

---

## Abstraction Is Not Always Win

Abstraction is a knife.  
In skilled hands, it cuts away noise.  
In clumsy hands, it obscures truth.

**Do not abstract prematurely.**  
Wait until you have three examples, not one.  
Wait until the pattern screams at you.  
Wait until the duplication *hurts*.

The wrong abstraction is worse than no abstraction.  
It is a lie that lives in your codebase,  
growing tendrils, demanding worship.

> *"Duplication is far cheaper than the wrong abstraction."*  
> — Sandi Metz

If your abstraction has more than two parameters,  
it is not an abstraction — it is a configuration file pretending to be code.

**Concrete is honest. Abstract is a promise. Keep promises sparingly.**

---

## Make Code Be Comments Instead of Making Comments Be Code

Comments are apologies.  
They are confessions that your code did not speak clearly enough.

**Good code whispers its intent in every line.**

```python
# Bad: comments explaining what
# Increment counter by 1
counter = counter + 1

# Good: code explaining why
users_who_completed_onboarding += 1
```

Name your variables like you're writing a letter to your future self.  
Name your functions like you're answering a question.

```
# A question
is_valid_email(string)

# An action
send_notification_to_user(user)

# A transformation
convert_celsius_to_fahrenheit(temp)
```

If you must comment, comment *why*, never *what*.  
The *what* should be obvious.  
If it's not, rewrite the code.

> *The best comment is the one you deleted because you fixed the code.*

---

## Never Nest

Nesting is cognitive tax.  
Every level of indentation is a wall between the reader and understanding.

**Flatten ruthlessly.**

```python
# Hell
if user:
    if user.is_active:
        if user.has_permission:
            if not user.is_banned:
                do_thing()

# Heaven
if not user:
    return
if not user.is_active:
    return
if not user.has_permission:
    return
if user.is_banned:
    return
    
do_thing()
```

**Guard clauses are elegance.**  
Fail fast. Fail early. Fail at the boundary.

> *"Arrow code" points right.  
> Good code flows down.*

If you need more than 3 levels of indentation, you have failed.  
Refactor. Extract. Invert the logic.

---

## Best Mode — Dependency Injection

**Hard-coded dependencies are hard-coded decisions.**  
They are concrete where you need flexibility.  
They are walls where you need doors.

Pass dependencies in.  
Do not reach out and grab them.

```python
# Coupled: the function owns its dependencies
def process_data():
    db = Database()
    api = ThirdPartyAPI()
    ...

# Decoupled: the function receives what it needs
def process_data(db, api):
    ...
```

**This is not about frameworks.**  
This is about *freedom*.

- You can test without touching the network.  
- You can swap implementations without rewriting logic.  
- You can compose without collision.

> *Make dependencies explicit.  
> Make implicit things impossible.*

If your function secretly depends on global state,  
it is not a function — it is a trap.

---

## Dear Functional Bros... I See You, I Hear You, and I Love You

**Functions are truth.**

A function is a contract:  
*Given this input, I promise this output.*  
No secrets. No surprises. No side effects hiding in closets.

### Pure Functions Are Sacred

```haskell
-- Same input, same output, always
f(x) = x * 2

-- No exceptions. No mutations. No lies.
```

**Purity is debuggability.**  
If a function is pure, you can test it in isolation.  
You can reason about it without context.  
You can trust it.

### Immutability Is Sanity

Mutation is action at a distance.  
You change something here, and something breaks over there.

**Do not mutate. Transform.**

```javascript
// Mutation: dangerous
let user = { name: "Alice", age: 30 }
user.age = 31  // Who else is watching this object?

// Transformation: safe
const user = { name: "Alice", age: 30 }
const olderUser = { ...user, age: 31 }
```

### Composition Is the Heartbeat

```
f ∘ g ∘ h

Small functions. Clear boundaries. Infinite combinations.
```

This is not dogma.  
This is *clarity crystallized into syntax*.

> *Object-oriented programming makes code understandable by encapsulating moving parts.  
> Functional programming makes code understandable by eliminating moving parts.*

---

## The Things We Must Add

### Names Are Spells

**A name summons meaning from the void.**

Bad names are curses that haunt your codebase.  
Good names are incantations that clarify intent.

```python
# Weak
def proc(d):
    return d * 1.2

# Strong  
def calculate_price_with_tax(price_before_tax):
    TAX_RATE = 1.2
    return price_before_tax * TAX_RATE
```

**Rules:**
- Booleans ask questions: `is_ready`, `has_permission`, `can_edit`
- Functions are verbs: `calculate`, `transform`, `validate`
- Classes are nouns: `User`, `EmailValidator`, `PaymentProcessor`
- Constants scream: `MAX_RETRY_ATTEMPTS`, `API_TIMEOUT_SECONDS`

One-letter variables are permitted in three cases only:
1. Loop indices: `i`, `j`, `k`
2. Mathematical conventions: `x`, `y`, `n`
3. Throwaway lambdas: `x => x * 2`

**Everywhere else: use real names.**

---

### Errors Are Not Exceptions

**Errors are data.**  
They should be handled like data, not thrown like grenades.

```python
# Exceptions: control flow in disguise
try:
    result = do_risky_thing()
except SomeError:
    # now what?

# Results: explicit success or failure
result = do_risky_thing()
if result.is_error():
    handle_error(result.error)
else:
    use_value(result.value)
```

Make failure *visible*.  
Make failure *expected*.  
Make failure *part of the type*.

```rust
// Rust knows this
fn divide(a: f64, b: f64) -> Result<f64, String> {
    if b == 0.0 {
        Err("Division by zero".to_string())
    } else {
        Ok(a / b)
    }
}
```

> *Exceptions are goto in disguise.*

---

### Timing Is a Side Effect

**If your code depends on *when* it runs, it is not a function — it is a spell bound to the moon.**

```python
# Time-dependent: fragile
def get_discount():
    if datetime.now().hour < 12:
        return 0.2
    return 0.1

# Time-agnostic: testable
def get_discount(current_hour):
    if current_hour < 12:
        return 0.2
    return 0.1
```

Pass time in.  
Do not reach for it.

**The same applies to:**
- Random numbers (pass the seed)
- Environment variables (pass the config)
- Global state (there is no global state)

---

## The Final Law

**Code is read 10 times for every time it is written.**

Optimize for the reader.  
Optimize for the future human staring at your work at 3 AM,  
exhausted, confused, desperate.

That human is you.

> *Write code as if the person who will maintain it is a violent psychopath who knows where you live.*

---

## Closing Invocation

There is no perfect code.  
There is only code that clearly expresses intent,  
handles failure gracefully,  
and can be understood by tired humans.

**Strive for:**
- Clarity over cleverness  
- Simplicity over sophistication  
- Explicitness over magic  

And when in doubt, **delete**.

---

*End of Manifesto*
